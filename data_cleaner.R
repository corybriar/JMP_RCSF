
pacman::p_load(tidyverse, ipumsr, maps, tidycensus, sf, tigris, tmap, tigris,
               raster, spData, usethis, stargazer, stringr, latex2exp, rvest, jsonlite, 
               httr, googleComputeEngineR, lfe, huxtable, stargazer, pracma, systemfit)
setwd("/Users/bigfoot13770/Documents/UO ECON PROGRAM/FUTILITY/JMP_RCSF")

options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

census_api_key("76eed5e2f67081d04c13b12ae20b0695817b7d68",
               install = TRUE, overwrite = TRUE)

# -------------IPUMS DATA ------------------

# Read in SF Shape
metro <- core_based_statistical_areas(cb=T, year = 2015) %>%
  filter(grepl("San Francisco", NAME))

# obtain puma shapes
pumas <- st_read("ipums_puma_2010") %>% st_transform(st_crs(metro))
pumas <- pumas[st_intersects(pumas,metro, sparse = F),] 


# Read in IPUMS data, filter down to Seattle
ddi <- read_ipums_ddi("usa_00012.xml")
sf_ipums <- read_ipums_micro(ddi) 
sf_ipums <- sf_ipums %>% filter(MET2013 == 41860, AGE > 18, !grepl("M", INDNAICS), INDNAICS != 0) %>%
  mutate(INDNAICS = str_sub(INDNAICS,1,2),
         skill = ifelse(EDUCD>100,"H","L"),
         OWNERSHP = OWNERSHP-1) 


# ---------TIDYCENSUS DATA---------------

acsvar <- load_variables(2019,"acs5")
# Read in spatial education data from the census
sanfran19 <- get_acs(geography = "tract",
                   variables = c(
                     pop25 = "B15002_001",
                     bachelorsM = "B15002_015",
                     bachelorsF = "B15002_032",
                     mastersM = "B15002_016",
                     mastersF = "B15002_033",
                     professionalM = "B15002_017", 
                     professionalF = "B15002_034", 
                     phdM = "B15002_018",
                     phdF = "B15002_035"),
                   year = 2019,
                   state = "CA",
                   output = "wide",
                   geometry = T) %>% 
  mutate(Hskill19 = bachelorsME +bachelorsFE +mastersME +mastersFE + professionalME + professionalFE + phdME + phdFE) %>%
  mutate(Lskill19 = pop25E - Hskill19) %>%
  dplyr::select(GEOID, NAME, Hskill19, Lskill19, geometry)

sanfran10 <- get_acs(geography = "tract",
                   variables = c(
                     pop25 = "B15002_001",
                     bachelorsM = "B15002_015",
                     bachelorsF = "B15002_032",
                     mastersM = "B15002_016",
                     mastersF = "B15002_033",
                     professionalM = "B15002_017", 
                     professionalF = "B15002_034", 
                     phdM = "B15002_018",
                     phdF = "B15002_035"),
                   year = 2010,
                   state = "CA",
                   output = "wide",
                   geometry = T) %>% 
  mutate(Hskill10 = bachelorsME +bachelorsFE +mastersME +mastersFE + professionalME + professionalFE + phdME + phdFE) %>%
  mutate(Lskill10 = pop25E - Hskill10) %>%
  dplyr::select(GEOID, NAME, Hskill10, Lskill10, geometry) %>%
  mutate(geometry = st_centroid(geometry))
# Merge dfs for 2019 and 2010
sanfran <- st_join(sanfran19,sanfran10) %>%
  dplyr::select(GEOID.x, NAME.x, Hskill19, Lskill19, Hskill10, Lskill10, geometry) %>%
  rename(GEOID = GEOID.x, NAME = NAME.x)
# Filter down to metro area
sanfran <- sanfran[st_intersects(sanfran$geometry,metro,sparse = F), ]

sf_shape <- get_acs(geography = "county",
                     variables = c(
                       pop25 = "B15002_001"),
                     year = 2010,
                     state = "CA",
                     output = "wide",
                     geometry = T) %>%
  filter(NAME == "San Francisco County, California")
sf_city <- sanfran[st_intersects(st_centroid(sanfran$geometry), sf_shape, sparse = F),]

#-------------BUILDING PERMIT DATA----------------------
# Permits from open data
res <- c("apartments", "1 family dwelling", "2 family dwelling")
build <- read_csv("Building_Permits.csv") %>% 
  mutate(`Street Name` = ifelse(substr(`Street Name`,1,1) == "0", substr(`Street Name`,2,4),`Street Name`)) %>%
  mutate(address = paste(`Street Number`,`Street Name`,sep = " "),
         cost = ifelse(is.na(`Revised Cost`) == T, `Estimated Cost`, `Revised Cost`),
         Location = gsub("POINT \\(", "XXX", Location)) %>%
  mutate(Location = gsub(" ", "YYY", Location)) %>%
  filter(`Existing Units` != `Proposed Units`,
         as.numeric(substr(`Current Status Date`,7,11)) >2009,
         as.numeric(substr(`Current Status Date`,7,11)) < 2020,
         `Existing Use` %in% res | `Proposed Use` %in% res) %>%
  rename(geometry = Location) %>%
  mutate(lon = as.numeric(gsub(".*XXX(.+)YYY.*", "\\1", geometry)),
         lat = as.numeric(gsub(".*YYY(.+)\\)", "\\1", geometry))) %>%
  filter(is.na(lon) == F, is.na(lat) == F) %>% st_as_sf(coords = c("lon","lat"))
  
st_set_crs(build, st_crs(sf_city)) 

# Write sf object for vm-rc to crank of for a couple days
st_write(build, "build_match_tract.shp")
build <- st_read("build_tract.shp")


#-----------PARCEL DATA-------------#
# Read in from VM
sf_city$parcels <- read_csv("parcel_tract.csv")$value


#-----------PUMA SORT---------------
pumas <- pumas[st_intersects(st_centroid(pumas),sf_shape,sparse = F), ]

puma_tract <- mapply(function(t){
  p <- pumas[st_intersects(st_centroid(t),pumas,sparse=F),]$PUMA[1]
  print(p)
  return(p)
},sf_city$geometry)

sf_city$puma <- puma_tract
sf_city <- sf_city %>% filter(Hskill19 >0)



























#------ Permit Files from Assessor's office

# Read in permit data
permitfiles <- dir(path = "permits/", pattern = ".xls")
# Filter out May2010 cause of a typo in STATUS_DATE
permitfiles <- c(permitfiles[1:72],permitfiles[74:108])

# Read in permit data
permits <- map_df(permitfiles, function(x){
  print(x)
  df <- read_excel(paste0("permits/",x)) 
  if("STREET_NUMBER" %in% colnames(df)){
    df <- df %>% dplyr::select(`APPLICATION #`, STATUS_DATE, `ESTIMATED COST`, `REVISED COST`, 
                                                       `EXISTING USE`, `EXISTING UNITS`, `PROPOSED USE`, `PROPOSED UNITS`,
                                                       STREET_NUMBER, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK)%>%
    mutate(year = parse_number(x), LOT = as.character(LOT))}
  else{
    df <- df %>% dplyr::select(`APPLICATION #`, STATUS_DATE, `ESTIMATED COST`, `REVISED COST`, 
                                                       `EXISTING USE`, `EXISTING UNITS`, `PROPOSED USE`, `PROPOSED UNITS`,
                                                       STREET_NUMBER..23, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>%
    mutate(year = parse_number(x), LOT = as.character(LOT))}
  return (df)
}) %>% 
  mutate(STREET_NUMBER = ifelse(year > 2013, STREET_NUMBER..23, STREET_NUMBER)) %>% 
  dplyr::select(-STREET_NUMBER..23)
# Prepare May2010 separately
may2010 <- read_excel(paste0("permits/May2010Issued_permits.xls")) %>% 
  dplyr::select(`APPLICATION #`, STATUS_DATE = `STATUS_ DATE`, `ESTIMATED COST`, `REVISED COST`, 
                `EXISTING USE`, `EXISTING UNITS`, `PROPOSED USE`, `PROPOSED UNITS`,
                STREET_NUMBER = `STREET_ NUMBER`, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% mutate(year = 2010, LOT = as.character(LOT))
# Combine with permits
permits <- rbind(permits,may2010)

res <- c("1 FAMILY DWELLING", "2 FAMILY DWELLING", "APARTMENTS")
# Filter down to relevant conversions
permits <- permits %>% rename(app_number = `APPLICATION #`, status_date = STATUS_DATE, est_cost = `ESTIMATED COST`, rev_cost = `REVISED COST`, 
                              ex_use = `EXISTING USE`, ex_units =`EXISTING UNITS`, prop_use = `PROPOSED USE`, 
                              prop_units = `PROPOSED UNITS`, street_number = STREET_NUMBER, street_name = AVS_STREET_NAME, 
                              street_sfx = AVS_STREET_SFX, description = DESCRIPTION, lot = LOT, block = BLOCK) %>%
  filter(prop_units != ex_units,
         prop_use %in% res | ex_use %in% res,
         est_cost != 1, rev_cost != 1,
         grepl("#",description) == F, grepl(" PA ", description) == F,
         grepl(" REF ", description) == F,
         is.na(street_number) == F) %>%
  mutate(deck = ifelse(grepl("deck", description) == T | grepl("DECK", description) == T | grepl("Deck", description) == T, 1,0),
         window = ifelse(grepl("window", description) == T | grepl("WINDOW", description) == T | grepl("Window", description) == T, 1,0),
         roof = ifelse(grepl("roof", description) == T | grepl("ROOF", description) == T | grepl("Roof", description) == T, 1,0),
         sprinkler = ifelse(grepl("sprinkler", description) == T | grepl("SPRINKLER", description) == T | grepl("Sprinkler", description) == T, 1,0),
         merge = ifelse(grepl("merge", description) == T | grepl("MERGE", description) == T | grepl("Merge", description) == T, 1,0),
         add = ifelse(grepl("add", description) == T | grepl("ADD", description) == T | grepl("Add", description) == T, 1,0),
         accessory = ifelse(grepl("accessory", description) == T | grepl("ACCESSORY", description) == T | grepl("Accessory", description) == T |
                              grepl("auxillary", description) == T | grepl("AUXILLARY", description) == T | grepl("Aukillary", description) == T |
                            grepl("ADU", description) == T, 1,0),
         remove = ifelse(grepl("remove", description) == T | grepl("REMOVE", description) == T | grepl("Remove", description) == T, 1,0),
         horizontal = ifelse(grepl("horizontal", description) == T | grepl("HORIZONTAL", description) == T | grepl("Horizontal", description) == T, 1,0),
         vertical = ifelse(grepl("vertical", description) == T | grepl("VERTICAL", description) == T | grepl("Vertical", description) == T, 1,0),
         unit_cost = rev_cost/abs(ex_units-prop_units)) %>%
  filter(unit_cost > 1000) %>%
  mutate(street_name = ifelse(substr(street_name,1,1) == "0",gsub("0","",street_name), street_name)) %>%
  mutate(address = paste(street_number,street_name, street_sfx, sep = " "),
         blocklot = paste(block, lot, sep = "-"))
  
write_csv(permits, "permits.csv")

# Find demolition files
demofiles <- dir(path = "demos/", pattern = ".xls")
demofiles <- c(demofiles[1:80], demofiles[82:101], demofiles[103:120])

demos <- map_df(demofiles, function(x){
  print(x)
  df <- read_excel(paste0("demos/",x)) 
  if("STREET_NUMBER" %in% colnames(df)){
    df <- df %>% dplyr::select(`APPLICATION #`, STATUS_DATE, `EXISTING USE`, `EXISTING UNITS`, `ESTIMATED COST`, `REVISED COST`, 
                  STREET_NUMBER, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% 
      mutate(year = parse_number(x), LOT = as.character(LOT), BLOCK = as.character(BLOCK))}
  else{
    df <- df %>% dplyr::select(`APPLICATION #`, STATUS_DATE, `EXISTING USE`, `EXISTING UNITS`, `ESTIMATED COST`, `REVISED COST`,
                               STREET_NUMBER..23, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% 
      mutate(year = parse_number(x), LOT = as.character(LOT), BLOCK = as.character(BLOCK))}
  return(df)
}) %>% mutate(STREET_NUMBER = ifelse(is.na(STREET_NUMBER..23), STREET_NUMBER, STREET_NUMBER..23)) %>%
  dplyr::select(-STREET_NUMBER..23)
# Change names for outliers 
may2010demo <- read_excel("demos/May2010Demolition_Permits.xls") %>%
  dplyr::select(`APPLICATION #`, STATUS_DATE, `EXISTING USE`, `EXISTING UNITS`, `ESTIMATED COST`, `REVISED COST`,
                STREET_NUMBER = `STREET_ NUMBER`, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% 
  mutate(year = 2010, LOT = as.character(LOT), BLOCK = as.character(BLOCK))
oct2010demo <- read_excel("demos/October2010Demolition_Permits.xls") %>%
  dplyr::select(`APPLICATION #`, STATUS_DATE, `EXISTING USE`, `EXISTING UNITS`, `ESTIMATED COST`, `REVISED COST`,
                STREET_NUMBER = `STREET_ NUMBER`, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% 
  mutate(year = 2010, LOT = as.character(LOT),BLOCK = as.character(BLOCK))
# Merge
demos <- rbind(demos, may2010demo) %>% rbind(oct2010demo)

#Filter down to relevant observations
demos <- demos %>% rename(app_number = `APPLICATION #`, status_date = STATUS_DATE, 
                          ex_use = `EXISTING USE`, ex_units = `EXISTING UNITS`,
                          street_number = STREET_NUMBER, street_name = AVS_STREET_NAME, 
                          street_sfx = AVS_STREET_SFX, description = DESCRIPTION, lot = LOT, block = BLOCK) %>%
  #filter(ex_use %in% res) %>%
  mutate(street_name = ifelse(substr(street_name,1,1) == "0",gsub("0","",street_name), street_name)) %>%
  mutate(address = paste(street_number,street_name, street_sfx, sep = " "),
         blocklot = paste(block, lot, sep = "-"))
# Save for viewing in Excel
write_csv(demos, "demos.csv") 

permits <- permits %>% mutate(conversion = ifelse(address %in% demos$address | blocklot %in% demos$blocklot, 1, 0))
  
  
  
  
  
  
  
  
  
  


