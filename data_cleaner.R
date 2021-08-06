
pacman::p_load(tidyverse, ipumsr, maps, tidycensus, sf, tigris, tmap, tigris,
               raster, spData, usethis, stargazer, stringr, latex2exp, rvest, jsonlite, 
               httr, googleComputeEngineR, lfe, huxtable, stargazer, pracma, systemfit, janitor, readxl)
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
ddi <- read_ipums_ddi("usa_00013.xml")
sf_ipums <- read_ipums_micro(ddi) 
sf_ipums <- sf_ipums %>% filter(MET2013 == 41860, AGE > 18, !grepl("M", INDNAICS), INDNAICS != 0, YEAR %in% c(2010,2019)) %>%
  mutate(INDNAICS = str_sub(INDNAICS,1,2),
         skill = ifelse(EDUCD>100,"H","L"),
         OWNERSHP = OWNERSHP-1) 
# Compute employment by industry 
sf_ipums <- sf_ipums %>% filter(CITY == 6290) %>% group_by(skill,YEAR) %>%
  mutate(skill_pop = n()) %>% 
  group_by(skill,INDNAICS) %>%
  mutate(ind_share = n()/skill_pop) %>% ungroup()
# Create df just of skill type, industry, and employment share
sf_indshares <- sf_ipums %>% filter(YEAR ==2019)
  distinct(skill,INDNAICS, .keep_all =T) %>%
  dplyr::select(skill, INDNAICS, ind_share)
# Determine occupancy by PUMA
sf_ipums <- sf_ipums %>%
  mutate(eta = ifelse(OWNERSHP == 1, "O",
                      ifelse(OWNERSHP == 0 & UNITSSTR > 3 & BUILTYR2 <6,"C", "X"))) %>%
  mutate(eta5ago = ifelse(eta == "C" &  MOVEDIN > 3, "C", "Not C"))

# filter by year, compute relevant moments
ipums10 <- sf_ipums %>% filter(YEAR == 2010) %>%
  group_by(PUMA) %>%
  mutate(CrateH10 = ifelse(eta == "C" & eta5ago == "Not C" & skill == "H",1,0),
         CrateL10 = ifelse(eta == "C" & eta5ago == "Not C" & skill == "L",1,0),
         CrateH10_C = ifelse(eta == "C" & eta5ago == "C" & skill == "H",1,0),
         CrateL10_C = ifelse(eta == "C" & eta5ago == "C" & skill == "L",1,0),
         XrateH10 = ifelse(eta == "X" & skill == "H",1,0),
         XrateL10 = ifelse(eta == "X" & skill == "L",1,0),
         OrateH10 = ifelse(eta == "O" & skill == "H",1,0),
         OrateL10 = ifelse(eta == "O" & skill == "L",1,0)
         ) %>%
  ungroup() %>% group_by(PUMA,skill) %>%
  mutate(CrateH10 = sum(CrateH10)/n(),
         CrateL10 = sum(CrateL10)/n(),
         CrateH10_C = sum(CrateH10_C)/n(),
         CrateL10_C = sum(CrateL10_C)/n(),
         XrateH10 = sum(XrateH10)/n(),
         XrateL10 = sum(XrateL10)/n(),
         OrateH10 = sum(OrateH10)/n(),
         OrateL10 = sum(OrateL10)/n()
  ) %>% distinct(PUMA, skill, .keep_all = T) %>% ungroup() %>%
  group_by(PUMA) %>%
  mutate(CrateH10 = sum(CrateH10),
                CrateL10 = sum(CrateL10),
                CrateH10_C = sum(CrateH10_C),
                CrateL10_C = sum(CrateL10_C),
                XrateH10 = sum(XrateH10),
                XrateL10 = sum(XrateL10),
                OrateH10 = sum(OrateH10),
                OrateL10 = sum(OrateL10)
  ) %>% distinct(PUMA, .keep_all = T) 

ipums19 <- sf_ipums %>% filter(YEAR == 2019) %>%
  group_by(PUMA) %>%
  mutate(CrateH19 = ifelse(eta == "C" & eta5ago == "Not C" & skill == "H",1,0),
         CrateL19 = ifelse(eta == "C" & eta5ago == "Not C" & skill == "L",1,0),
         CrateH19_C = ifelse(eta == "C" & eta5ago == "C" & skill == "H",1,0),
         CrateL19_C = ifelse(eta == "C" & eta5ago == "C" & skill == "L",1,0),
         XrateH19 = ifelse(eta == "X" & skill == "H",1,0),
         XrateL19 = ifelse(eta == "X" & skill == "L",1,0),
         OrateH19 = ifelse(eta == "O" & skill == "H",1,0),
         OrateL19 = ifelse(eta == "O" & skill == "L",1,0)
  ) %>%
  ungroup() %>% group_by(PUMA,skill) %>%
  mutate(CrateH19 = sum(CrateH19)/n(),
         CrateL19 = sum(CrateL19)/n(),
         CrateH19_C = sum(CrateH19_C)/n(),
         CrateL19_C = sum(CrateL19_C)/n(),
         XrateH19 = sum(XrateH19)/n(),
         XrateL19 = sum(XrateL19)/n(),
         OrateH19 = sum(OrateH19)/n(),
         OrateL19 = sum(OrateL19)/n()
  ) %>% distinct(PUMA, skill, .keep_all = T) %>% ungroup() %>%
  group_by(PUMA) %>%
  mutate(CrateH19 = sum(CrateH19),
         CrateL19 = sum(CrateL19),
         CrateH19_C = sum(CrateH19_C),
         CrateL19_C = sum(CrateL19_C),
         XrateH19 = sum(XrateH19),
         XrateL19 = sum(XrateL19),
         OrateH19 = sum(OrateH19),
         OrateL19 = sum(OrateL19)
  ) %>% distinct(PUMA, .keep_all = T)


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
  filter(pop25E > 0) %>%
  mutate(Lskill19 = pop25E - Hskill19) %>%
  dplyr::select(GEOID, NAME, Hskill19, Lskill19, geometry) 

sanfran19

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
  filter(pop25E >0) %>%
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

#-----------PUMA SORT---------------
pumas <- pumas[st_intersects(st_centroid(pumas),sf_shape,sparse = F), ]

puma_tract <- mapply(function(t){
  p <- pumas[st_intersects(st_centroid(t),pumas,sparse=F),]$PUMA[1]
  print(p)
  return(p)
},sf_city$geometry)

sf_city$puma <- puma_tract
sf_city$puma[is.na(sf_city$puma) == T] <- "07502"

# Merge with Ipums data
puma_merge <- map_df(sf_city$NAME,function(x){
  tract <- sf_city %>% filter(NAME == x)
  p <- tract$puma[1] %>% as.character() %>% substr(5,5)
  p10 <- ipums10[substr(as.character(ipums10$PUMA),4,4) == p,]
  p19 <- ipums19[substr(as.character(ipums19$PUMA),4,4) == p,]
  df <- tibble(
    NAME = tract$NAME[1],
    CrateH10 = p10$CrateH10[1],
    CrateL10 = p10$CrateL10[1],
    CrateH10_C = p10$CrateH10_C[1],
    CrateL10_C = p10$CrateL10_C[1],
    XrateH10 = p10$XrateH10[1],
    XrateL10 = p10$XrateL10[1],
    OrateH10 = p10$OrateH10[1],
    OrateL10 = p10$OrateL10[1],
    CrateH19 = p19$CrateH19[1],
    CrateL19 = p19$CrateL19[1],
    CrateH19_C = p19$CrateH19_C[1],
    CrateL19_C = p19$CrateL19_C[1],
    XrateH19 = p19$XrateH19[1],
    XrateL19 = p19$XrateL19[1],
    OrateH19 = p19$OrateH19[1],
    OrateL19 = p19$OrateL19[1]
  )
})

# Merge into sf_city
sf_city <- sf_city %>% left_join(puma_merge, by = "NAME") %>% distinct(NAME, .keep_all = T) 

sf_city <- sf_city %>% 
  mutate(cent = st_centroid(geometry),
         CshareH10 = CrateH10*Hskill10/sum(sanfran$Hskill10, na.rm = T),
         CshareL10 = CrateL10*Lskill10/sum(sanfran$Lskill10, na.rm = T),
         CshareH10_C = CrateH10_C*Hskill10/sum(sanfran$Hskill10, na.rm = T),
         CshareL10_C = CrateL10_C*Lskill10/sum(sanfran$Lskill10, na.rm = T),
         XshareH10 = XrateH10*Hskill10/sum(sanfran$Hskill10, na.rm = T),
         XshareL10 = XrateL10*Lskill10/sum(sanfran$Lskill10, na.rm = T),
         OshareH10 = OrateH10*Hskill10/sum(sanfran$Hskill10, na.rm = T),
         OshareL10 = OrateL10*Lskill10/sum(sanfran$Lskill10, na.rm = T),
         
         CshareH19 = CrateH19*Hskill19/sum(sanfran$Hskill19),
         CshareL19 = CrateL19*Lskill19/sum(sanfran$Lskill19),
         CshareH19_C = CrateH19_C*Hskill19/sum(sanfran$Hskill19),
         CshareL19_C = CrateL19_C*Lskill19/sum(sanfran$Lskill19),
         XshareH19 = XrateH19*Hskill19/sum(sanfran$Hskill19),
         XshareL19 = XrateL19*Lskill19/sum(sanfran$Lskill19),
         OshareH19 = OrateH19*Hskill19/sum(sanfran$Hskill19),
         OshareL19 = OrateL19*Lskill19/sum(sanfran$Lskill19)
         )

#-------------BUILDING PERMIT DATA----------------------
# Permits from open data
res <- c("apartments", "1 family dwelling", "2 family dwelling")
build <- read_csv("Building_Permits.csv") %>% 
  mutate(`Street Name` = ifelse(substr(`Street Name`,1,1) == "0", substr(`Street Name`,2,4),`Street Name`)) %>%
  mutate(address = paste(`Street Number`,`Street Name`,sep = " "),
         cost = ifelse(is.na(`Revised Cost`) == T, `Estimated Cost`, `Revised Cost`),
         Location = gsub("POINT \\(", "XXX", Location)) %>%
  mutate(Location = gsub(" ", "YYY", Location)) %>%
  filter(as.numeric(substr(`Current Status Date`,7,11)) >2009,
         as.numeric(substr(`Current Status Date`,7,11)) < 2020,
         `Existing Use` %in% res | `Proposed Use` %in% res) %>%
  rename(geometry = Location) %>%
  mutate(lon = as.numeric(gsub(".*XXX(.+)YYY.*", "\\1", geometry)),
         lat = as.numeric(gsub(".*YYY(.+)\\)", "\\1", geometry))) %>%
  filter(is.na(lon) == F, is.na(lat) == F,`Current Status` == "complete") %>% st_as_sf(coords = c("lon","lat"))
                                                        
build_ref <- mapply(function(x){
  num <- regmatches(x, gregexpr("[[:digit:]]+", x)) %>% 
    unlist() %>%  as.numeric() 
  ref <- sum(ifelse(num > 1990 & num < 2020,1,0))
  
}, build$Description)

build <- build %>% mutate(ref = as.character(build_ref)) %>%
  filter(ref == 0, grepl("revision", Description) == F) %>%
  mutate(merge = ifelse(grepl("merge",Description) == T| grepl("combine",Description) == T, 1, 0),
         sprinklers = ifelse(grepl("sprinkler", Description) == T | grepl("ssystem", Description) == T, 1, 0),
         vert = ifelse(grepl("vertical",Description) == T, 1, 0),
         conversion = ifelse(grepl("conversion",Description) == T, 1, 0),
         office = ifelse(grepl("office",Description) == T, 1, 0),
         legal = ifelse(grepl("legalize",Description) == T, 1, 0),
         horizontal = ifelse(grepl("horizontal",Description) == T, 1, 0),
         parking = ifelse(grepl("parking",Description) == T, 1, 0),
         condo = ifelse(grepl("condo",Description) == T, 1, 0),
         address = toupper(paste(address, `Street Suffix`, " ")),
         `Existing Units` = ifelse(is.na(`Existing Use`) == T, 0, `Existing Units`),
         blocklot = paste(gsub("^0+", "", gsub("^0+", "", Block)), gsub("^0+", "", gsub("^0+", "", Lot)), sep = "-"))
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
                               STREET_NUMBER...23, AVS_STREET_NAME, AVS_STREET_SFX, DESCRIPTION, LOT, BLOCK) %>% 
      mutate(year = parse_number(x), LOT = as.character(LOT), BLOCK = as.character(BLOCK))}
  return(df)
}) %>% mutate(STREET_NUMBER = ifelse(is.na(STREET_NUMBER...23), STREET_NUMBER, STREET_NUMBER...23)) %>%
  dplyr::select(-STREET_NUMBER...23)
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

demos <- read_csv("demos.csv") %>% mutate(cost = ifelse(is.na(`REVISED COST`), `ESTIMATED COST`, `REVISED COST`)) 



# Cross reference with demolition data 
demolitions <- map_df(1:dim(build)[1],function(x){
  dfa <- demos %>% filter(address == build$address[x])
  dfb <- demos %>% filter(blocklot == build$blocklt[x])
  if(dim(dfa)[1] == 0 & dim(dfb)[1] == 0){
    df <- tibble(
      demo_units = 0,
      demo_costs = 0
    )
  } else if(dim(dfa)[1] == 0){
    df <- tibble(
      demo_units = sum(dfb$ex_units, na.rm = T),
      demo_costs = sum(dfb$cost, na.rm = T)
    )
  } else {
    df <- tibble(
      demo_units = sum(dfa$ex_units, na.rm = T),
      demo_costs = sum(dfa$cost, na.rm = T)
    )
  }
  return(df)
})


build <- build %>% 
  filter(ExstngUn != PrpsdUn | condo ==1 | convrsn == 1| demo >0 | PrmtTyD =="demolitions") %>%
  mutate(demo = ifelse(demo == 2, 1, demo)) %>%
  mutate(demo = ifelse(PrmtTyD == "demolition",1,demo)) %>%
  cbind(demolitions)

st_set_crs(build, st_crs(sf_city)) 

# Write sf object for vm-rc to crank of for a couple days
st_write(build, "build_match_tract.shp")

build <- st_read("build_match_tract.shp")

# Read back in
build <- st_read("build_tract.shp")


# Read in data on condo conversions and evictions
evictions <- st_read("Eviction Notices/evictions.shp") %>% st_set_crs(st_crs(sf_city))

conlotto <- read_csv("condo_lottery.csv", col_types = "ccccdddc") %>%
  mutate(blocklot = paste(gsub("^0+", "", gsub("^0+", "", Block)), gsub("^0+", "", gsub("^0+", "", Lot)), sep = "-"))

conwin <- conlotto %>% filter(Status == "WINNER")
# Compare lotto winners to build dataframe
lotterywin <- mapply(function(x){
  win <- ifelse(x %in% conwin$blocklot,1,0)
},build$blocklt)

build$lotto <- lotterywin

build <- st_set_crs(build, st_crs(sf_city))

# Manipulate conlotto to get conditional probabilities of winning
conlotto <- conlotto %>% 
  mutate(win = ifelse(Status == "WINNER",1,0)) %>%
  group_by(Number_Of_Units) %>%
  mutate(prob = sum(win)/n()) %>%
  ungroup() %>% group_by(Lottery_Year, Status) %>%
  mutate(struct_num = n(),
         units_convo = n()*Number_Of_Units) %>% ungroup() 

zip_shape <- get_acs(geography = "zip code tabulation area",
                     variables = c(
                       pop25 = "B15002_001"),
                     year = 2019,
                     output = "wide",
                     geometry = T) 
zip_shape <- zip_shape[st_intersects(st_centroid(zip_shape),sf_shape,sparse = F),]

# Assign tracts to zip codes by majority share
zip_tract <- map_df(1:195, function(x){
  tract <- sf_city[x,]
  zip <- zip_shape[st_intersects(zip_shape,tract, sparse = F),][which.max(st_area(st_intersection(zip_shape,tract))),]
  df <- tibble(
    NAME = tract$NAME[1],
    main_zip = zip$GEOID[1]
  )
  return(df)
})

sf_city <- sf_city %>% left_join(zip_tract, by = "NAME")

tract_share <- map_df(unique(sf_city$main_zip), function(z){
  fake_zip <- sf_city %>% filter(main_zip == z)
  df <- tibble(
    NAME = fake_zip$NAME,
    zip_share = as.numeric(st_area(fake_zip$geometry)/st_area(st_union(fake_zip)))
  )
  return(df)
})

sf_city <- sf_city %>% left_join(tract_share, by = "NAME")

conwin <- conlotto %>% filter(win ==1) %>%
  group_by(Zip_Code,Number_Of_Units) %>%
  mutate(winners = n()) %>%
  distinct(Zip_Code,Number_Of_Units, .keep_all = T) %>%
  mutate(Zip_Code = substr(Zip_Code,1,5)) %>%
  filter(Zip_Code %in% zip_shape$GEOID)

conwin_zip <- map_df(unique(zip_shape$GEOID), function(z){
  zip <- filter(conwin, Zip_Code == z)
  df <- tibble(
    main_zip = z,
    win2=zip[zip$Number_Of_Units == 2,]$winners[1],
    win3=zip[zip$Number_Of_Units == 3,]$winners[1],
    win4=zip[zip$Number_Of_Units == 4,]$winners[1],
    win5=zip[zip$Number_Of_Units == 5,]$winners[1],
    win6=zip[zip$Number_Of_Units == 6,]$winners[1],
  )
})

conwin_zip[is.na(conwin_zip)] <- 0

sf_city <- sf_city %>% left_join(conwin_zip, by = "main_zip") %>%
  mutate(win2 = win2 * zip_share/sum(conwin_zip$win2),
         win3 = win3 * zip_share/sum(conwin_zip$win3),
         win4 = win4 * zip_share/sum(conwin_zip$win4),
         win5 = win5 * zip_share/sum(conwin_zip$win5),
         win6 = win6 * zip_share/sum(conwin_zip$win6))

#------Land Use Data-------
landuse <- st_read("Land Use/landuse.shp") %>% st_set_crs(st_crs(build)) %>%
  mutate(blocklot = paste(gsub("^0+", "", gsub("^0+", "", block_num)), gsub("^0+", "", gsub("^0+", "", lot_num)), sep = "-"))
st_write(build, "build_year.shp")

st_write(landuse, "landuse_year.shp")

# yrbuilt <- mapply(function(x){
#   yearbuilt <- landuse[st_intersects(x,landuse$geometry, sparse = F),]$yrbuilt[1]
#   
# }, build$geometry[1:10])
yrbuilt <- read_csv("yearbuilt.csv")
build$yrbuilt <- yrbuilt$value




#-----------PARCEL DATA-------------#
# Read in from VM
sf_city$parcels <- read_csv("parcel_tract.csv")$value



sf_city <- sf_city %>% filter(Hskill19 >0)




#------------NHGIS DATA--------------

#------Units Data----------

units19 <- read_csv("NHGIS_data/housing19.csv") %>% filter(NAME_E %in% sf_city$NAME) %>%
  mutate(Cstock19 = AMJXE069 + AMJXE070 +AMJXE071 +AMJXE072 +AMJXE076 +AMJXE077 +AMJXE078 +AMJXE079 +
           AMJXE083 +AMJXE084 +AMJXE085 +AMJXE086,
         Cstruct19 = round(AMJXE069/3 + AMJXE070/12 +AMJXE071/34.5 +AMJXE072/50 +AMJXE076/3 +AMJXE077/12 +AMJXE078/34.5 +AMJXE079/50 +
           AMJXE083/3 +AMJXE084/12 +AMJXE085/34.5 +AMJXE086/50,0),
         Xstock19 = AMJXE047 +AMJXE048 +AMJXE049 +AMJXE050 +AMJXE051 +AMJXE054 +AMJXE055 +AMJXE056 + AMJXE057 +AMJXE058 +
           AMJXE061 + AMJXE062 + AMJXE063 + AMJXE064 + AMJXE065 + AMJXE068 +AMJXE075 +AMJXE082,
         Xstruct19 = round(AMJXE047 +AMJXE048/3 +AMJXE049/12 +AMJXE050/34.5 +AMJXE051/50 +AMJXE054 +AMJXE055/3 +AMJXE056/12 + AMJXE057/34.5 +AMJXE058/50 +
           AMJXE061 + AMJXE062/3 + AMJXE063/12 + AMJXE064/34.5 + AMJXE065/50 + AMJXE068 +AMJXE075 +AMJXE082,0),
         Ostock19 = AMJXE002,
         Ostruct19 = round(AMJXE004 + AMJXE005/3 + AMJXE006/12 + AMJXE007/34.5 + AMJXE008/50 + 
           AMJXE011 + AMJXE012/3 + AMJXE013/12 + AMJXE014/34.5 + AMJXE015/50 + 
           AMJXE018 + AMJXE019/3 + AMJXE020/12 + AMJXE021/34.5 + AMJXE022/50 + 
           AMJXE025 + AMJXE026/3 + AMJXE027/12 + AMJXE028/34.5 + AMJXE029/50 + 
           AMJXE032 + AMJXE033/3 + AMJXE034/12 + AMJXE035/34.5 + AMJXE036/50 + 
           AMJXE039 + AMJXE040/3+ AMJXE041/12 + AMJXE042/34.5 + AMJXE043, 0)) %>%
  rename(NAME = NAME_E) 

nanner <- function(x){
  ifelse(is.nan(x) == T, 0, x)
}


units10 <- read_csv("NHGIS_data/housing10.csv") %>% filter(NAME_E %in% sf_city$NAME)%>%
  mutate(Cstock10 = J9GE055 + J9GE056 + J9GE057 + J9GE058 + J9GE062 + J9GE063 + J9GE064 + J9GE065 + 
           J9GE069 + J9GE070 + J9GE071 + J9GE072,
         Cstruct10 = round(J9GE055/3 + J9GE056/14.5 + J9GE057/34.5 + J9GE058/50 + J9GE062/3 + J9GE063/14.5 + J9GE064/34.5 + J9GE065/50 + 
           J9GE069/3 + J9GE070/14.5 + J9GE071/34.5 + J9GE072/50,0),
         Xstock10 = J9GE040 + J9GE041 + J9GE042 + J9GE043 + J9GE044 + J9GE047 + J9GE048 + J9GE049 + J9GE050 + J9GE051 +
           J9GE054 + J9GE061 + J9GE068,
         Xstruct10 = round(J9GE040 + J9GE041/3 + J9GE042/14.5 + J9GE043/34.5 + J9GE044/50 + J9GE047 + J9GE048/3 + J9GE049/14.5 + J9GE050/34.5 + J9GE051/50 +
           J9GE054 + J9GE061 + J9GE068,0),
         Ostock10 = J9GE002,
         Ostruct10 = round(J9GE004 + J9GE005/3 + J9GE006/12 + J9GE007/34.5 + J9GE008/50 +
           J9GE011 + J9GE012/3 + J9GE013/12 + J9GE014/34.5 + J9GE015/50 +
           J9GE018 + J9GE019/3 + J9GE020/12 + J9GE021/34.5 + J9GE022/50 +
           J9GE025 + J9GE026/3 + J9GE027/12 + J9GE028/34.5 + J9GE039/50 +
           J9GE032 + J9GE033/3 + J9GE034/12 + J9GE035/34.5  + J9GE036/50,0)) %>%
  mutate(Csize10 = Cstock10/Cstruct10, Xsize10=Xstock10/Xstruct10, Osize10 = Ostock10/Ostruct10) %>%
  mutate(Csigma10 = sqrt((3-Csize10)^2*J9GE055/3 + (12-Csize10)^2*J9GE056/14.5 + (34.5-Csize10)^2*J9GE057/34.5 + (50-Csize10)^2*J9GE058/50 + 
                         (3-Csize10)^2*J9GE062/3 + (12-Csize10)^2*J9GE063/14.5 + (34.5-Csize10)^2*J9GE064/34.5 + (50-Csize10)^2*J9GE065/50 + 
                         (3-Csize10)^2*J9GE069/3 + (12-Csize10)^2*J9GE070/14.5 + (34.5-Csize10)^2*J9GE071/34.5 + (50-Csize10)^2*J9GE072/50),
         Xsigma10 = sqrt((1-Xsize10)^2*J9GE040 + (3-Xsize10)^2*J9GE041/3 + (12-Xsize10)^2*J9GE042/14.5 + (34.5-Xsize10)^2*J9GE043/34.5 + (50-Xsize10)^2*J9GE044/50 + 
                         (1-Xsize10)^2*J9GE047 + (3-Xsize10)^2*J9GE048/3 + (12-Xsize10)^2*J9GE049/14.5 + (34.5-Xsize10)^2*J9GE050/34.5 + (50-Xsize10)^2*J9GE051/50 +
                         (1-Xsize10)^2*J9GE054 + (1-Xsize10)^2*J9GE061 + (1-Xsize10)^2*J9GE068),
         Osigma10 = sqrt((1-Osize10^2)*J9GE004 + (3-Osize10)^2*J9GE005 + (12-Osize10)^2*J9GE006 + (34.5-Osize10)*J9GE007 + (50-Osize10)^2*J9GE008 +
           (1-Osize10^2)*J9GE011 + (3-Osize10)^2*J9GE012 + (12-Osize10)^2*J9GE013 + (34.5-Osize10)*J9GE014 + (50-Osize10)^2*J9GE015 +
           (1-Osize10^2)*J9GE018 + (3-Osize10)^2*J9GE019 + (12-Osize10)^2*J9GE020 + (34.5-Osize10)*J9GE021 + (50-Osize10)^2*J9GE022 +
           (1-Osize10^2)*J9GE025 + (3-Osize10)^2^2*J9GE026 + (34.5-Osize10)*J9GE027 + (50-Osize10)^2*J9GE029 +
           (1-Osize10^2)*J9GE032 + (3-Osize10)^2*J9GE033 + (12-Osize10)^2*J9GE034 + (34.5-Osize10)*J9GE035 + (50-Osize10)^2*J9GE036)) %>%
  rename(NAME = NAME_E) 
        
#------ Rent data ---------
na_means <- function(x,y){
  ifelse(is.na(mean(c(x,y), na.rm = T)) == T, 0, mean(c(x,y), na.rm = T))
}

unitrents10 <- read_csv("NHGIS_data/rents/rents10.csv") %>% filter(NAME %in% sf_city$NAME) %>% 
  right_join(units10, by = "NAME") %>%
  mutate(
    Crent10 = (na_means(b1970,b1960)*(J9GE055 + J9GE056 + J9GE057 + J9GE058) + 
      na_means(b1950,b1940)*(J9GE062 + J9GE063 + J9GE064 + J9GE065) + 
      na_means(b1939,b1939)*(J9GE069 + J9GE070 + J9GE071 + J9GE072))/Cstock10,
    Xrent10 = (na_means(b2005,b2000)*(J9GE040 + J9GE041 + J9GE042 + J9GE043 + J9GE044) + 
      na_means(b1990,b1980)*(J9GE047 + J9GE048 + J9GE049 + J9GE050 +J9GE051)+ 
      na_means(b1970,b1960)*J9GE054 + na_means(b1950,b1940)*J9GE061 + na_means(b1939,b1939)*J9GE068)/Xstock10
  ) %>%
  mutate(
    CsigmaR10 = sqrt(nanner((na_means(b1970_me,b1960_me)*sqrt(J9GE055 + J9GE056 + J9GE057 + J9GE058)/qt(0.975,J9GE055 + J9GE056 + J9GE057 + J9GE058))^2) +
      nanner((na_means(b1950_me,b1940_me)*sqrt(J9GE062 + J9GE063 + J9GE064 + J9GE065)/qt(0.975,J9GE062 + J9GE063 + J9GE064 + J9GE065))^2) +
      nanner((na_means(b1939_me,b1939_me)*sqrt(J9GE069 + J9GE070 + J9GE071 + J9GE072)/qt(0.975,J9GE069 + J9GE070 + J9GE071 + J9GE072))^2)), 
    XsigmaR10 = sqrt(nanner((na_means(b2005_me,b2000_me)*sqrt(J9GE040 + J9GE041 + J9GE042 + J9GE043 + J9GE044)/qt(0.975,J9GE040 + J9GE041 + J9GE042 + J9GE043 + J9GE044))^2) +
      nanner((na_means(b1990_me,b1980_me)*sqrt(J9GE047 + J9GE048 + J9GE049 + J9GE050 + J9GE051)/qt(0.975,J9GE047 + J9GE048 + J9GE049 + J9GE050 + J9GE051))^2) +
      nanner((na_means(b1970_me,b1960_me)*sqrt(J9GE054)/qt(0.975,J9GE054))^2) + nanner((na_means(b1950_me,b1940_me)*sqrt(J9GE061)/qt(0.975,J9GE061))^2) + 
      nanner((na_means(b1939_me,b1939_me)*sqrt(J9GE068)/qt(0.975,J9GE068))^2))
  )

unitrents19 <- read_csv("NHGIS_data/rents/rents19.csv") %>% filter(NAME %in% sf_city$NAME) %>% 
  right_join(units19, by = "NAME") %>%
  mutate(
    Crent19 = (na_means(b1970,b1960)*(AMJXE069 + AMJXE070 +AMJXE071 +AMJXE072) + 
                 na_means(b1950,b1940)*(AMJXE076 + AMJXE077 +AMJXE078 +AMJXE079) + 
                 na_means(b1939,b1939)*(AMJXE083 + AMJXE084 +AMJXE085 +AMJXE086))/Cstock19,
    Xrent19 = (na_means(b2014,b2010)*(AMJXE047 + AMJXE048 +AMJXE049 +AMJXE050 + AMJXE051) + 
                 na_means(b2000,b2000)*(AMJXE054 + AMJXE055 +AMJXE056 +AMJXE057 + AMJXE058) + 
                 na_means(b1990,b1980)*(AMJXE061 + AMJXE062 +AMJXE063 +AMJXE064 + AMJXE065) +
                 na_means(b1970,b1960)*AMJXE068 + na_means(b1950,b1940)*AMJXE075 + na_means(b1939,b1939)*AMJXE082)/Xstock19
  )%>%
  mutate(
    CsigmaR19 = sqrt(nanner((na_means(b1970_me,b1960_me)*sqrt(AMJXE069 + AMJXE070 +AMJXE071 +AMJXE072)/qt(0.975,AMJXE069 + AMJXE070 +AMJXE071 +AMJXE072))^2) +
                       nanner((na_means(b1950_me,b1940_me)*sqrt(AMJXE076 + AMJXE077 +AMJXE078 +AMJXE079)/qt(0.975,AMJXE076 + AMJXE077 +AMJXE078 +AMJXE079))^2) +
                       nanner((na_means(b1939_me,b1939_me)*sqrt(AMJXE083 + AMJXE084 +AMJXE085 +AMJXE086)/qt(0.975,AMJXE083 + AMJXE084 +AMJXE085 +AMJXE086))^2)), 
    XsigmaR19 = sqrt(nanner((na_means(b2014_me,b2000_me)*sqrt(AMJXE047 + AMJXE048 +AMJXE049 +AMJXE050 + AMJXE051)/qt(0.975,AMJXE047 + AMJXE048 +AMJXE049 +AMJXE050 + AMJXE051))^2) +
                       nanner((na_means(b1990_me,b1980_me)*sqrt(AMJXE061 + AMJXE062 +AMJXE063 +AMJXE064 + AMJXE065)/qt(0.975,AMJXE061 + AMJXE062 +AMJXE063 +AMJXE064 + AMJXE065))^2) +
                       nanner((na_means(b1970_me,b1960_me)*sqrt(AMJXE068)/qt(0.975,AMJXE068))^2) + nanner((na_means(b1950_me,b1940_me)*sqrt(AMJXE075)/qt(0.975,AMJXE075))^2) + 
                       nanner((na_means(b1939_me,b1939_me)*sqrt(AMJXE082)/qt(0.975,AMJXE082))^2))
  )

# Pair down dataframes to merge with unitrents
merge19 <- unitrents19 %>% dplyr::select(NAME, Cstock19, Cstruct19, Crent19, Xstock19, Xstruct19, Xrent19, Ostock19, Ostruct19)
merge10 <- unitrents10 %>% 
  dplyr::select(NAME, Cstock10, Csize10, Cstruct10, Crent10, Csigma10, CsigmaR10,Xstock10, Xsize10, Xstruct10, Xrent10, Xsigma10, XsigmaR10, Ostock10, Osize10, Ostruct10, Osigma10)

sf_city <- left_join(sf_city, merge19, by = "NAME") %>% left_join(merge10,by = "NAME")


#--------------Mortgage data -----------------#
mortgages19 <- read_csv("NHGIS_data/mortgages/mortgages19.csv") %>% 
  filter(NAME_E %in% sf_city$NAME) %>% rename(NAME = NAME_E) %>%
  mutate(mortgage19 = (200*AL1QE003 + 250*AL1QE004 + 350*AL1QE005 +450*AL1QE006 +550*AL1QE007 +
           650*AL1QE008 +750*AL1QE009 +850*AL1QE010 +950*AL1QE011 +1125*AL1QE012 +1375*AL1QE013 +
           1750*AL1QE014 +2250*AL1QE015 +2750*AL1QE016 +3250*AL1QE017 +3750*AL1QE018 + 4000*AL1QE019)/
           (AL1QE003 + AL1QE004 + AL1QE005 + AL1QE006 + AL1QE007 + AL1QE008 + AL1QE009 + AL1QE010 + AL1QE011 +
              AL1QE012 + AL1QE013 + AL1QE014 + AL1QE015 + AL1QE016 + AL1QE017 + AL1QE018 + AL1QE019)) %>%
  mutate(Mortgagesigma19 = sqrt((mortgage19-200)^2*AL1QE003 + (mortgage19-250)^2*AL1QE004 + 
           (mortgage19-350)^2*AL1QE005 +(mortgage19-450)^2*AL1QE006 +(mortgage19-550)^2*AL1QE007 +
           (mortgage19-650)^2*AL1QE008 +(mortgage19-750)^2*AL1QE009 +(mortgage19-850)^2*AL1QE010 +
           (mortgage19-950)^2*AL1QE011 +(mortgage19-1125)^2*AL1QE012 +(mortgage19-1375)^2*AL1QE013 +
             (mortgage19-1750)^2*AL1QE014 +(mortgage19-2250)^2*AL1QE015 +(mortgage19-2750)^2*AL1QE016 +
             (mortgage19-3250)^2*AL1QE017 +(mortgage19-3750)^2*AL1QE018 + (mortgage19-4000)^2*AL1QE019)) %>%
  dplyr::select(NAME,mortgage19,Mortgagesigma19)

mortgages10 <- read_csv("NHGIS_data/mortgages/mortgages10.csv") %>%
  filter(NAME %in% sf_city$NAME) %>% 
  mutate(mortgage10 = (200*JTRE003 + 250*JTRE004 + 350*JTRE005 +450*JTRE006 +550*JTRE007 +
                         650*JTRE008 +750*JTRE009 +850*JTRE010 +950*JTRE011 +1125*JTRE012 +1375*JTRE013 +
                         1750*JTRE014 +2250*JTRE015 +2750*JTRE016 +3250*JTRE017 +3750*JTRE018 + 4000*JTRE019)/
           (JTRE003 + JTRE004 + JTRE005 +JTRE006 +JTRE007 +JTRE008 +JTRE009 +JTRE010 +JTRE011 +JTRE012 +JTRE013 +
              JTRE014 + JTRE015 + JTRE016 + JTRE017 + JTRE018 + JTRE019)) %>%
  mutate(Mortgagesigma10 = sqrt((mortgage10-200)^2*JTRE003 + (mortgage10-250)^2*JTRE004 + 
                                  (mortgage10-350)^2*JTRE005 +(mortgage10-450)^2*JTRE006 +(mortgage10-550)^2*JTRE007 +
                                  (mortgage10-650)^2*JTRE008 +(mortgage10-750)^2*JTRE009 +(mortgage10-850)^2*JTRE010 +
                                  (mortgage10-950)^2*JTRE011 +(mortgage10-1125)^2*JTRE012 +(mortgage10-1375)^2*JTRE013 +
                                  (mortgage10-1750)^2*JTRE014 +(mortgage10-2250)^2*JTRE015 +(mortgage10-2750)^2*JTRE016 +
                                  (mortgage10-3250)^2*JTRE017 +(mortgage10-3750)^2*JTRE018 + (mortgage10-4000)^2*JTRE019)) %>%
  dplyr::select(NAME,mortgage10,Mortgagesigma10)

sf_city <- sf_city %>% left_join(mortgages10, by = "NAME") %>% left_join(mortgages19, by = "NAME")

# buyout data
buyouts <- st_read("Buyout Agreements/buyouts.shp") %>%
  filter(is.na(number_of_) == F, is.na(buyout_amo) == F) %>% st_set_crs(st_crs(sf_city)) %>%
  filter(as.numeric(substr(as.character(date_buyou),1,4)) < 2020 & 
           as.numeric(substr(as.character(date_buyou),1,4)) >2010)

buyout_tracts <- mapply(function(x){
  temp <- sf_city[st_intersects(sf_city$geometry,x, sparse = F),]$NAME[1]
}, buyouts$geometry) 

buyouts <- buyouts %>% 
  mutate(NAME = buyout_tracts, buy_unit = buyout_amo/number_of_) %>% group_by(NAME) %>%
  mutate(avgBO = mean(buy_unit)) %>% ungroup() %>% distinct(NAME, .keep_all = T) %>%
  dplyr::select(NAME,avgBO) 

buyouts$geometry <- NULL

sf_city <- left_join(sf_city,buyouts, by = "NAME") 
sf_city <- sf_city %>%  mutate(avgBO = ifelse(is.na(avgBO) == T, mean(avgBO,na.rm = T),avgBO))

#--------------TRACT CHARACTERISTICS---------------

# Central business district
cbd <- st_sfc(st_point(c(-122.4029,37.7952)), crs = st_crs(sf_city))

# Read in ZBP data from 2018
zbp18 <- read_csv("zbp18detail.txt") %>% 
  filter(grepl("----",naics), !grepl("------",naics)) %>%
  mutate(`n<5` = as.numeric(ifelse(`n<5` == "N", 0, `n<5`)),
         n5_9 = as.numeric(ifelse(n5_9 == "N", 0, n5_9)),
         n10_19 = as.numeric(ifelse(n10_19 == "N", 0, n10_19)),
         n20_49 = as.numeric(ifelse(n20_49 == "N", 0, n20_49)),
         n50_99 = as.numeric(ifelse(n50_99 == "N", 0, n50_99)),
         n100_249 = as.numeric(ifelse(n100_249 == "N", 0, n100_249)),
         n250_499 = as.numeric(ifelse(n250_499 == "N", 0, n250_499)),
         n500_999 = as.numeric(ifelse(n500_999 == "N", 0, n500_999)),
         n1000 = as.numeric(ifelse(n1000 == "N", 0, n1000)),
         naics = substr(naics, 1,2)) %>% 
  complete(zip, naics, fill = list(`n<5`=0,n5_9=0,n10_19=0,n20_49=0,n50_99=0,n100_249=0,n250_499=0,n500_999=0,n1000=0,est=0)) %>%
  mutate(total = 3*`n<5` + 7*n5_9 + 14.5*n10_19 + 34.5*n20_49 + 74.5*n50_99 + 174.5*n100_249 + 374.5*n250_499 + 749.5*n500_999 + 1000*n1000)
  
# Pull zipcodes for SF, use to filter down data frame before "rotating"
zip_shape <- get_acs(geography = "zip code tabulation area",
                    variables = c(
                      pop25 = "B15002_001"),
                    year = 2019,
                    output = "wide",
                    geometry = T) 
zip_shape <- zip_shape[st_intersects(st_centroid(zip_shape),sf_shape,sparse = F),]
zbp18 <- zbp18 %>% filter(zip %in% zip_shape$GEOID)

# Create variable for each sector, "rotate" data frame into zbp
zbp <- map_df(substr(as.character(unique(zip_shape$GEOID)), 1, 5),function(z){
  temp <- zbp18 %>% filter(zip == z)
  df = st_sf(
    zip = z,
    total11 = temp[temp$naics == "11",]$total[1],
    total21 = temp[temp$naics == "21",]$total[1],
    total22 = temp[temp$naics == "22",]$total[1],
    total23 = temp[temp$naics == "23",]$total[1],
    total31 = temp[temp$naics == "31",]$total[1],
    total42 = temp[temp$naics == "42",]$total[1],
    total44 = temp[temp$naics == "44",]$total[1],
    total48 = temp[temp$naics == "48",]$total[1],
    total51 = temp[temp$naics == "51",]$total[1],
    total52 = temp[temp$naics == "52",]$total[1],
    total53 = temp[temp$naics == "53",]$total[1],
    total54 = temp[temp$naics == "54",]$total[1],
    total55 = temp[temp$naics == "55",]$total[1],
    total56 = temp[temp$naics == "56",]$total[1],
    total61 = temp[temp$naics == "61",]$total[1],
    total62 = temp[temp$naics == "62",]$total[1],
    total71 = temp[temp$naics == "71",]$total[1],
    total72 = temp[temp$naics == "72",]$total[1],
    total81 = temp[temp$naics == "81",]$total[1],
    total99 = temp[temp$naics == "99",]$total[1],
    geometry = zip_shape[zip_shape$GEOID == as.numeric(z),]$geometry[1]
  )
}) %>% filter(is.na(total11) == F) %>%
  mutate(
    share11 = ifelse(is.na(total11/sum(total11))==T,0,total11/sum(total11)),
    share21 = ifelse(is.na(total21/sum(total21))==T,0,total21/sum(total21)),
    share22 = ifelse(is.na(total22/sum(total22))==T,0,total22/sum(total22)),
    share23 = ifelse(is.na(total23/sum(total23))==T,0,total23/sum(total23)),
    share31 = ifelse(is.na(total31/sum(total31))==T,0,total31/sum(total31)),
    share42 = ifelse(is.na(total42/sum(total42))==T,0,total42/sum(total42)),
    share44 = ifelse(is.na(total44/sum(total44))==T,0,total44/sum(total44)),
    share48 = ifelse(is.na(total48/sum(total48))==T,0,total48/sum(total48)),
    share51 = ifelse(is.na(total51/sum(total51))==T,0,total51/sum(total51)),
    share52 = ifelse(is.na(total52/sum(total52))==T,0,total52/sum(total52)),
    share53 = ifelse(is.na(total53/sum(total53))==T,0,total53/sum(total53)),
    share54 = ifelse(is.na(total54/sum(total54))==T,0,total54/sum(total54)),
    share55 = ifelse(is.na(total55/sum(total55))==T,0,total55/sum(total55)),
    share56 = ifelse(is.na(total56/sum(total56))==T,0,total56/sum(total56)),
    share61 = ifelse(is.na(total61/sum(total61))==T,0,total61/sum(total61)),
    share62 = ifelse(is.na(total62/sum(total62))==T,0,total62/sum(total62)),
    share71 = ifelse(is.na(total71/sum(total71))==T,0,total71/sum(total71)),
    share72 = ifelse(is.na(total72/sum(total72))==T,0,total72/sum(total72)),
    share81 = ifelse(is.na(total81/sum(total81))==T,0,total81/sum(total81)),
    share99 = ifelse(is.na(total99/sum(total99))==T,0,total99/sum(total99))
  ) %>% st_set_crs(st_crs(sf_city))

# Filter indshares to match zbp data
sf_indshares <- sf_indshares %>% filter(INDNAICS %in% as.character(unique(zbp18$naics))) %>%
  complete(skill,INDNAICS,fill = list(ind_share = 0))
# Extract Matricies to compute expected commute times
emp_shareH = as.matrix(sf_indshares$ind_share[1:20])
emp_shareL = as.matrix(sf_indshares$ind_share[21:40])
loc_shares = zbp[,23:42] %>% as.tibble() %>% dplyr::select(!geometry) %>% as.matrix()
dist <- st_distance(st_centroid(sf_city$geometry), st_centroid(zbp$geometry)) %>% as.numeric() %>% as.matrix() %>% Reshape(195, 27)

# Compute expected commute distances
sf_city$Hcommute <- dist %*% loc_shares %*% emp_shareH
sf_city$Lcommute <- dist %*% loc_shares %*% emp_shareL
sf_city$cbd <- st_distance(st_centroid(sf_city),cbd) %>% as.numeric()

rm(zip_shape, zbp18)

#---------Park Data----------
parks <- read_csv("nanda_parks_tract.csv") %>% filter(tract_fips10 %in% sf_city$GEOID) %>%
  dplyr::select(tract_fips10, count_open_parks, tot_park_area, prop_park_area_tract) %>%
  rename(GEOID = tract_fips10, parks = count_open_parks, park_area = tot_park_area, park_share = prop_park_area_tract)

sf_city <- left_join(sf_city, parks, by = "GEOID")

#--------Water Data---------
# Read in boundary file of CA, find distance of each tract to water
california <- get_acs(geography = "state",
                                  variables = c(
                                    pop25 = "B15002_001"),
                                  year = 2010,
                                  state = "CA",
                                  output = "wide",
                                  geometry = T) %>% st_boundary()
sf_city$water <- st_distance(st_centroid(sf_city$geometry),california) %>% as.numeric()
rm(california)

#------Zoning Data--------
sf_city <- sf_city %>% left_join(read_csv("zone_share.csv"), by = "NAME")






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



permits <- permits %>% mutate(conversion = ifelse(address %in% demos$address | blocklot %in% demos$blocklot, 1, 0))

# ----- REGRESSION ANALYSIS

# Filter build into adjustment/renovation data
adjustments <- build %>% mutate(ExstngUn = ifelse(is.na(ExstngUn) == T, 0, ExstngUn),
                                PrpsdUn = ifelse(is.na(PrpsdUn) == T,0,PrpsdUn)) %>%
  filter(merge == 0, PrmtTyp != "demolitions", demo == 0, ExstngUn >0, ExstngUs %in% res,
         ExstngUn < PrpsdUn, ExstngUn > 0, cost > 1) %>%
  mutate(unit_cost = cost/(PrpsdUn - ExstngUn),
         scaling = PrpsdUn/ExstngUn) %>%
  filter(unit_cost > 2000)

condocons <- build %>% 
  mutate(ExstngUn = ifelse(is.na(ExstngUn) == T, 0, ExstngUn),
         PrpsdUn = ifelse(is.na(PrpsdUn) == T,0,PrpsdUn),
         condocon = condo*convrsn) %>%
  filter(PrmtTyp != "demolitions", ExstngUn >0, ExstngUs %in% res,
         condocon == 1|lotto == 1, ExstngUn > 0, cost > 1) %>%
  mutate(unit_cost = cost/(PrpsdUn - ExstngUn)) %>%
  filter(unit_cost > 2000)

construction <- build %>% 
  mutate(ExstngUn = ifelse(is.na(ExstngUn) == T, 0, ExstngUn),
         PrpsdUn = ifelse(is.na(PrpsdUn) == T,0,PrpsdUn)) %>%
  filter(PrmtTyD != "demolitions", (demo ==1 | ExstngUn == 0) , PrpsdUn > ExstngUn) %>%
  mutate(unit_cost = cost/(PrpsdUn - ExstngUn)) %>% filter(unit_cost >2000)

app <- mapply(function(x){
  num <- sum(ifelse(as.numeric(substr(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))),1,4)) >2000 & 
                      as.numeric(substr(unlist(regmatches(x, gregexpr("[[:digit:]]+", x))),1,4)) <2020,1,0))
  return(num)
}, construction$Dscrptn)

construction <- construction %>% mutate(is_app = app) %>%
  filter(is_app == 0, grepl("new construction",PrmtTyD))

removals <- build %>% mutate(ExstngUn = ifelse(is.na(ExstngUn) == T, 0, ExstngUn),
                             PrpsdUn = ifelse(is.na(PrpsdUn) == T,0,PrpsdUn)) %>%
  filter(PrmtTyD != "demolitions", demo == 0, ExstngUn >0, ExstngUs %in% res,
         ExstngUn > PrpsdUn, ExstngUn > 0, cost > 1) %>%
  mutate(unit_cost = cost/(-PrpsdUn +ExstngUn)) %>%
  filter(unit_cost > 2000)

#-------ADJUSTMENT COSTS REGRESSION--------

regA <- lm(log(cost) ~ log(PrpsdUn/ExstngUn) + tract -1, adjustments)

FE_A <- tibble(
  GEOID = sort(unique(adjustments$tract)),
  Gamma = exp(regA$coefficients[2:167])*regA$coefficients[1]*(-2)/(-1)
)

sf_city <- left_join(sf_city,FE_A, by = "GEOID") %>%
  mutate(Gamma = ifelse(is.na(Gamma) == T, mean(Gamma, na.rm = T), Gamma))

regCondo <- lm(cost ~ PrpsdUn, condocons)

#-------RENOVATION REGRESSIONS----------

construction2 <- filter(construction, log(PrpsdUn) > 1)

regC <- lm(log(cost) ~ log(PrpsdUn) + tract -1, construction)

FE_C <- tibble(
  GEOID = sort(unique(construction$tract)),
  Psi = exp(regC$coefficients[2:149])*regC$coefficients[1]*(-2)/(-1)
)

sf_city <- left_join(sf_city, FE_C, by = "GEOID") %>% 
  mutate(Psi = ifelse(is.na(Psi) == T, mean(Psi, na.rm = T), Psi))

construction$rawfit<- exp(predict(regC))

removals <- mutate(removals, demo_units = ExstngUn-PrpsdUn)

regR <- lm(cost ~ demo_units + merge, removals)



#------- WRITE SF_CITY TO CSV -----





ddi <- read_ipums_ddi("usa_00013.xml")
ipums <- read_ipums_micro(ddi) 
ipums <- ipums %>% filter(MET2013 == 41860, AGE > 18, !grepl("M", INDNAICS), INDNAICS != 0, YEAR ==2019) %>%
  mutate(INDNAICS = str_sub(INDNAICS,1,2),
         skill = ifelse(EDUCD>100,"H","L"),
         OWNERSHP = OWNERSHP-1) 
  

neighborhoods <- st_read("neighborhoods/neighborhoods.shp") %>%
  rename(GEOID = geoid) %>%
  dplyr::select(GEOID,nhood)

neighborhoods$geometry <- NULL

winners <- conwin_zip %>% 
  mutate(win = win2 + win3 + win4 + win5 + win6) %>%
  dplyr::select(main_zip,win)

sf_city<- left_join(sf_city,winners, by = "main_zip") %>%
  group_by(main_zip) %>%
  mutate(win = as.numeric(win * st_area(geometry)/st_area(st_union(geometry)))) %>%
  ungroup() %>% mutate(win = win/sum(win))


sf_neighborhoods <- sf_city %>% left_join(neighborhoods, by = "GEOID") %>%
  group_by(nhood) %>% 
    mutate(
    Hskill19 = sum(Hskill19, na.rm = T),        
    Lskill19 = sum(Lskill19, na.rm = T),       
    Hskill10 = sum(Hskill10, na.rm = T),      
    Lskill10 = sum(Lskill10, na.rm = T),       
    CrateH10 = mean(CrateH10, na.rm = T),       
    CrateL10 = mean(CrateL10, na.rm = T),       
    CrateH10_C = mean(CrateH10_C, na.rm = T),    
    CrateL10_C = mean(CrateL10_C, na.rm = T),     
    XrateH10 = mean(XrateH10, na.rm = T),       
    XrateL10 = mean(XrateL10, na.rm = T),       
    OrateH10 = mean(OrateH10, na.rm = T),       
    OrateL10 = mean(OrateL10, na.rm = T),      
    CrateH19 = mean(CrateH19, na.rm = T),       
    CrateL19 = mean(CrateL19, na.rm = T),       
    CrateH19_C = mean(CrateH19_C, na.rm = T),     
    CrateL19_C = mean(CrateL19_C, na.rm = T),     
    XrateH19 = mean(XrateH19, na.rm = T),      
    XrateL19 = mean(XrateL19, na.rm = T),       
    OrateH19 = mean(OrateH19, na.rm = T),       
    OrateL19 = mean(OrateL19, na.rm = T),       
    CshareH10 = sum(CshareH10, na.rm = T),        
    CshareL10 = sum(CshareL10, na.rm = T),    
    CshareH10_C = sum(CshareH10_C, na.rm = T),    
    CshareL10_C = sum(CshareL10_C, na.rm = T),    
    XshareH10 = sum(XshareH10, na.rm = T),      
    XshareL10 = sum(XshareL10, na.rm = T),      
    OshareH10 = sum(OshareH10, na.rm = T),      
    OshareL10 = sum(OshareL10, na.rm = T),       
    CshareH19 = sum(CshareH19, na.rm = T),      
    CshareL19 = sum(CshareL19, na.rm = T),      
    CshareH19_C = sum(CshareH19_C, na.rm = T),    
    CshareL19_C = sum(CshareL19_C, na.rm = T),   
    XshareH19 = sum(XshareH19, na.rm = T),      
    XshareL19 = sum(XshareL19, na.rm = T),      
    OshareH19 = sum(OshareH19, na.rm = T),      
    OshareL19 = sum(OshareL19, na.rm = T),      
    win2 = sum(win2, na.rm = T),           
    win3 = sum(win3, na.rm = T),           
    win4 = sum(win4, na.rm = T),           
    win5 = sum(win5, na.rm = T),          
    win6 = sum(win6, na.rm = T), 
    win = sum(win),
    parcels = sum(parcels, na.rm = T),        
    Cstock19 = sum(Cstock19, na.rm = T),       
    Cstruct19 = sum(Cstruct19, na.rm = T),      
    Crent19  = mean(Crent19, na.rm = T),    
    Xstock19 = sum(Xstock19, na.rm = T),        
    Xstruct19 = sum(Xstruct19, na.rm = T),      
    Xrent19  = mean(Xrent19, na.rm = T),       
    Ostock19 = sum(Ostock19, na.rm = T),       
    Ostruct19 = sum(Ostruct19, na.rm = T),    
    Cstock10  = sum(Cstock10, na.rm = T),      
    Csize10  = mean(Csize10, na.rm = T),      
    Cstruct10 = sum(Cstruct10, na.rm = T),      
    Crent10  = mean(Crent10, na.rm = T),       
    Csigma10  = mean(Csigma10, na.rm = T),     
    CsigmaR10 = mean(CsigmaR10, na.rm = T),      
    Xstock10  = sum(Xstock10, na.rm = T),      
    Xsize10 = mean(Xsize10, na.rm = T),        
    Xstruct10  = sum(Xstruct10, na.rm = T),     
    Xrent10 = mean(Xrent10, na.rm = T),       
    Xsigma10 = mean(Xsigma10, na.rm = T),       
    XsigmaR10 = mean(XsigmaR10, na.rm = T),      
    Ostock10 = sum(Ostock10, na.rm = T),     
    Osize10 = mean(Osize10, na.rm = T),        
    Ostruct10 = sum(Ostruct10, na.rm = T),   
    Osigma10 = mean(Osigma10, na.rm = T),      
    mortgage10 = mean(mortgage10, na.rm = T),     
    Mortgagesigma10 = mean(Mortgagesigma10, na.rm = T),
    mortgage19 = mean(mortgage19, na.rm = T),     
    Mortgagesigma19 = mean(Mortgagesigma19, na.rm = T),
    avgBO = mean(avgBO, na.rm = T),          
    Hcommute = mean(Hcommute, na.rm = T),       
    Lcommute = mean(Lcommute, na.rm = T),       
    parks = sum(parks, na.rm = T),         
    park_area = sum(park_area, na.rm = T),      
    park_share = as.numeric(sum(park_area, na.rm = T)/st_area(st_union(geometry))),     
    water = as.numeric(st_distance(st_centroid(st_union(geometry)),california)),          
    Gamma = mean(Gamma, na.rm = T),          
    Psi = mean(Psi, na.rm = T), 
    geometry = st_union(geometry)
  ) %>%
  distinct(nhood, .keep_all = T) %>%
  dplyr::select(-cbd, -res_share, -com_share, -ind_share, -mixed_share, -park_share, -GEOID, -NAME)

# recalculate cbd, instruments for each neighbohood
leftovers <- map_df(1:dim(sf_neighborhoods)[1],function(x){
  hood <- sf_neighborhoods[x,]
  zoning <- zones[st_intersects(hood,zones, sparse = F),]
  c = st_set_crs(cbd,st_crs(hood))
  df = tibble(
    nhood = hood$nhood[1],
    park_share = as.numeric(hood$park_area/st_area(hood)),
    cbd = as.numeric(st_distance(c,st_centroid(hood))),
    res_share = ifelse(length(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Residential",]),0),hood))[1] == 0, 0,
      as.numeric(st_area(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Residential",]),0),hood)))/as.numeric(st_area(hood))),
    com_share = ifelse(length(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Commercial",]),0),hood))[1] == 0, 0,
        as.numeric(st_area(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Commercial",]),0),hood)))/as.numeric(st_area(hood))),
    ind_share = ifelse(length(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Industrial",]),0),hood))[1] ==0, 0,
      as.numeric(st_area(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Industrial",]),0),hood)))/as.numeric(st_area(hood))),
    mix_share = ifelse(length(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Mixed Use",]),0),hood))[1] ==0, 0,
      as.numeric(st_area(st_intersection(st_buffer(st_union(zoning[zoning$gen == "Mixed Use",]),0),hood)))/as.numeric(st_area(hood))),
  )
  return(df)
})

sf_neighborhoods <- sf_neighborhoods %>% left_join(leftovers, by = "nhood")

sf_data <- sf_neighborhoods
sf_data$geometry <- NULL
sf_data$cent <- NULL

sf_data[is.na(sf_data)] <- 0
write_csv(sf_data, "sf_data.csv")





