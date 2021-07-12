
pacman::p_load(tidyverse, ipumsr, maps, tidycensus, sf, tigris, tmap, tigris,
               raster, spData, usethis, stargazer, stringr, latex2exp, rvest, jsonlite, 
               httr, googleComputeEngineR, lfe, huxtable, stargazer, pracma, systemfit, janitor, fixedest)
setwd("/Users/bigfoot13770/Documents/UO ECON PROGRAM/FUTILITY/JMP_RCSF")

options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

census_api_key("76eed5e2f67081d04c13b12ae20b0695817b7d68",
               install = TRUE, overwrite = TRUE)

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

regCondo <- lm(cost ~ PrpsdUn, condocons)

#-------RENOVATION REGRESSIONS----------

construction2 <- filter(construction, log(PrpsdUn) > 2)

regC <- lm(log(cost) ~ log(PrpsdUn) + tract -1, construction)

construction$rawfit<- exp(predict(regC))

regR <- lm(cost ~ (ExstngUn-PrpsdUn) + merge, removals)




