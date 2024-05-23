# Magali Blanco
# 05/23/2024
# script to prepare files for analysis in 2_health_impact.R. output is 'analysis_files.rda' 

##################################################################################################
# NOTES 
##################################################################################################
# FILES NEEDED
# exposure (air pollution) data 
#    grid/area predictions
#    a study outline
# population data
# incidence/prevalence data
# health impact functions

#####################################
# resources
## nice visual of cumulative impacts approach: https://www.eea.europa.eu/publications/assessing-the-risks-to-health
## mortality data from CD WONDER database: http://wonder.cdc.gov 
## outcomes need to be combined w/ Census pop info?? (see manual Appendix D)
## Census data
### book: Analysing Census Data, by Walker: https://walker-data.com/census-r/index.html  
### presentation: https://walker-data.com/umich-workshop-2022/intro-2020-census/#14 

#####################################
# note that this originally was coded for BenMap. Some coding refernces this (e.g COL/ROW variables)

##################################################################################################
# SETUP
##################################################################################################
# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(data.table, # for read_aermod_excel()
               readxl, tidyverse, sf, 
               tigris, tidycensus,  # census data
               rstudioapi, #askForPassword()
               leaflet #map QC checks
)

source("functions.R")

# where modified data will be saved
dt_path <- file.path("data", "modified")
# where files for analysis will be saved
out_path <- file.path("data", "output", "analysis")

# create directories if they don't exist
directory_list <- c(dt_path,
                    out_path,
                    file.path("data", "raw", "cdc_wonder", "api_request"))

lapply(directory_list, function(x) {if(!dir.exists(x)) {dir.create(x, recursive = T)}})

set.seed(1)
##################################################################################################
# COMMON VARIABLES
##################################################################################################
# grid should be projected into NAD83 (or WGS84) if want to use in BenMAP easily
## project CRS:
lat_long_crs <- "+proj=longlat +datum=WGS84"  # same as: 4326 
utm_zone10_crs <- 32610 # crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +ellps=WGS84"
# m_crs <- 32148 # NAD83 Washington North

# census data year
yr <- 2020

# file path for predicted air pollution
ap_file_path <- file.path("~", "Documents", "AERMOD UFP Analysis", "WillThesisAERMODfiles", "AERMOD Import R")

##################################################################################################
# GRIDDED AIR POLLUTION PREDICTIONS & STUDY AREA
##################################################################################################
# function reads AERMOD excel file
read_aermod_excel <- function(filename, direction, activity){
  
  north_grid = lapply(1:7, FUN = function(x){
    north_grid = data.table(read_excel(filename, sheet= x))
    colnames(north_grid)[1] = "meters"
    rownames(north_grid) = north_grid$meters
    north_grid$meters= NULL
    
    new_grid= data.table(x = -999, y  = -999, value = -999)
    
    grid1 = rbindlist(lapply(1: nrow(north_grid),FUN = function(i) {
      rbindlist(lapply(1: ncol(north_grid), FUN = function(j) {
        #print(i)
        #print(j)
        # MB 3/29/24: changed x/y definitions, which were switched
        data.table(x= colnames(north_grid)[j],
                   y= rownames(north_grid)[i],
                   value = as.numeric(north_grid[i, ..j]))
      }))
    }))
    grid1
  })
  
  north_grid = rbindlist(north_grid)
  
  north_grid[, x := as.numeric(x)]
  north_grid[, y := as.numeric(y)]
  north_grid[, aermod := as.numeric(value)]
  north_grid$direction = direction
  north_grid$activity = activity
  north_grid$value = NULL
  return(north_grid)
}
##################################################################################################

# --> why does this grid have a small tail on bottom right of 9 centroids?

# file xlsx can be any aermod landing/take-off file w/ the same grid system
grid <- read_aermod_excel(file.path(ap_file_path, "NorthBoundLanding_Predictions.xlsx"),"NFlow","landing") %>%
  st_as_sf(coords = c('x', 'y'), crs=utm_zone10_crs, remove = F) %>% 
  cbind(st_coordinates(.)) %>%
  select(x=X, y=Y)  

# format grid based on the creating_a_regular_grid_for_use_with_benmap_0.docx from Henry Raab (BenMap contact)
# grid must have 'Column' (or Col) and 'Row'

# get the grid resolution. # 1000 m
grid_resolution <- grid %>% distinct(x) %>% arrange(x) %>% pull(x)
grid_resolution <- grid_resolution[2] - grid_resolution[1]
# # check that y is same. # looks good
# grid_resolution <- grid %>% distinct(y) %>% arrange(y) %>% pull(y)
# grid_resolution <- x_resolution[2] - x_resolution[1]

x_min <- min(grid$x)
y_min <- min(grid$y)

# amount that needs to be added so that each X/Y value becomes 1
col_adjustment <- 1- x_min/grid_resolution
row_adjustment <- 1- y_min/grid_resolution

# See Table 5.1 for columns necessary
grid <- grid %>%
  mutate(
    Column = x/grid_resolution + col_adjustment,
    Row = y/grid_resolution + row_adjustment) %>%
  st_transform(lat_long_crs) 

# save grid with original x/y & calibrate AERMOD predictions in aermod_calibration.R
saveRDS(grid, file.path(dt_path, "grid.rds"))

##################################################################################################
# AERMOD PREDICTIONS
# Import calibrated EARMOD predictions (from aermod_calibration.R)

# --> TO DO: upate this using Ningrui's work?
calibrated_aermod <- readRDS(file.path(ap_file_path, "data", "modified", "grid_aermod_predictions_calibrated.rda")) %>%
  st_drop_geometry()

grid_baseline <- left_join(grid, calibrated_aermod, by=c("Column", "Row")) %>%
  mutate(
    # using daily averages to estimate an annual average
    Metric = "D24HourMean",
    'Seasonal Metric' = "QuarterlyMean", #need to have this for e.g., like PM2.5?
    'Annual Metric' = "Mean") %>%
  select(Column, Row, Metric, `Seasonal Metric`, `Annual Metric`, Values) %>%
  mutate(long = st_coordinates(.)[1],
         lat = st_coordinates(.)[2]) #%>% st_drop_geometry()


##################################################################################################
# STUDY OUTLINE 
##################################################################################################
aermod_area <- grid %>%
  # make into multi-point 
  summarize() %>%
  # take the outline
  st_convex_hull() %>%
  mutate(Column = 1,
         Row = 1) 

#saveRDS(aermod_area, file.path(dt_path, "aermod_area.rda"))

##################################################################################################
# CENSUS TRACTS (DEMOGRAPHIC & AGGREGATION UNITS)
##################################################################################################
study_counties <- c("King County", "Snohomish County", "Pierce County", "Kitsap County", "Island County")

tract_path <- file.path("data", "raw", "wa_kc_census_tracts.rda")
if(file.exists(tract_path)) {
  tracts0 <- readRDS(tract_path)
} else {
  tracts0 <- tracts(state = "WA", county = study_counties, year=yr) %>%
    filter(ALAND!=0) %>%
    st_transform(crs = lat_long_crs) 
  saveRDS(tracts0, tract_path)
}

##################################################################################################
# function returns areas within the study area
areas_overlap <- function(x., y.=aermod_area) {
 
  intersects <- x. %>% 
    #st_transform(crs = utm_zone10_crs) %>%
    st_intersects(x =., y=y.) %>% 
    as.data.frame() %>%
    pull(row.id)
  
  x.[intersects,] %>% #plot()
    select(GEOID) %>%
    # unique identifier (for BenMAP results to be linked)
    mutate(Column = 1, 
           Row = row_number())
}


##################################################################################################
# in study area
tracts <- areas_overlap(x.=tracts0, y.=aermod_area)  

# leaflet() %>%
#   addTiles() %>%
#   addPolygons(data=blocks) %>%
#   addPolygons(data=study_area)

# ggplot() + 
#   geom_sf(data=study_area, aes(fill="study")) + 
#   geom_sf(data=blocks, aes(fill=POP20), alpha=0.3)  

##################################################################################################
# AIR POLLUTION AT AGGREGATION AREAS (E.G., CENSUS TRACTS)
##################################################################################################
# #QC: visualize data
# ggplot() + 
#   geom_sf(data=tracts) + 
#   geom_sf(data=grid_baseline, aes(col=Values), size=1, alpha=0.5)

# average air pollution predictions if there are multple within an aggregation area
multiple_predictions_per_ag_location <- grid_baseline %>%
  select(-c(Row, Column, long, lat)) %>%
  st_join(tracts, .) %>%
  group_by(GEOID) %>%
  mutate(ap_points_per_location =n (),
         ap_points_per_location = ifelse(is.na(Values), 0, ap_points_per_location)) %>% 
  group_by(across(c(-Values))) %>%
  summarize(Values=mean(Values)) %>%
  ungroup()

# for locations without AP values, take the nearest estimate
no_ap_location <- multiple_predictions_per_ag_location %>%
  filter(is.na(Values)) 

nearest_ap_location <- no_ap_location %>% 
  st_nearest_feature(grid_baseline) %>%
  st_drop_geometry()

nearest_value <- grid_baseline[nearest_ap_location, "Values"] %>% st_drop_geometry() %>% pull()

no_ap_location <- no_ap_location %>%
  mutate(Values = nearest_value)  

air_pollution <- multiple_predictions_per_ag_location %>%
  drop_na(Values) %>%
  bind_rows(no_ap_location)

# QC: check that things look right. looks good!
# ggplot(data=air_pollution) + geom_sf(aes(fill=Values))

##################################################################################################
# CENSUS POPULATION DATA
##################################################################################################
# decennial variables are at the TRACT level; ACS are at the BLOCK GROUP level & averaged for 5 yr periods
# For the decennial Census, possible dataset choices include "pl" for the redistricting files; "dhc" for the Demographic and Housing Characteristics file and "dp" for the Demographic Profile (2020 only), and "sf1" or "sf2" (2000 and 2010) and "sf3" or "sf4" (2000 only) for the various summary files. 

# api_key <- askForPassword("Enter your Census API Key.\nTo create one, see: https://api.census.gov/data/key_signup.html") 
# census_api_key(api_key, install=T)

# --> do we trust the 2020 census? 

## demographics
decennial_variables_dp <- load_variables(year = yr, dataset = "dp") %>%  
  filter(#grepl("Count", label),
    label %in% c("Count!!SEX AND AGE!!Total population",
                 "Count!!SEX AND AGE!!Male population",
                 "Count!!SEX AND AGE!!Female population",
                 "Count!!MEDIAN AGE BY SEX!!Both sexes",
                 "Count!!SEX AND AGE!!Total population!!Selected Age Categories!!18 years and over",
                 "Count!!SEX AND AGE!!Total population!!Selected Age Categories!!65 years and over") | 
      (grepl("Count!!SEX AND AGE!!Total population!!", label) & grepl("years$", label)) |
      label == "Count!!RACE!!Total population" | grepl("Count!!RACE!!Total population!!One Race!!", label) |  label == "Count!!RACE!!Total population!!Two or More Races" |
      grepl("Count!!HISPANIC OR LATINO!!Total population!!", label)  | #probably won't actually use "white" w/o ethnicity?
      label %in% c(#"Count!!HISPANIC OR LATINO BY RACE!!Total population",
        "Count!!HISPANIC OR LATINO BY RACE!!Total population!!Not Hispanic or Latino!!White alone" #NH White
      ) |
      grepl("Count!!RELATIONSHIP!!Total population!!In households!!", label) & grepl("sex", label)) %>%
  mutate(dataset = "dp",
         year = yr,
         variable = case_when(
           name == "DP1_0001C" ~ "pop_total",
           
           name == "DP1_0002C" ~ "age_00_04",
           name == "DP1_0003C" ~ "age_05_09",
           name == "DP1_0004C" ~ "age_10_14",
           name == "DP1_0005C" ~ "age_15_19",
           name == "DP1_0006C" ~ "age_20_24",
           name == "DP1_0007C" ~ "age_25_29",
           name == "DP1_0008C" ~ "age_30_34",
           name == "DP1_0009C" ~ "age_35_39",
           name == "DP1_0010C" ~ "age_40_44",
           name == "DP1_0011C" ~ "age_45_49",
           name == "DP1_0012C" ~ "age_50_54",
           name == "DP1_0013C" ~ "age_55_59",
           name == "DP1_0014C" ~ "age_60_64",
           name == "DP1_0015C" ~ "age_65_69",
           name == "DP1_0016C" ~ "age_70_74",
           name == "DP1_0017C" ~ "age_75_79",
           name == "DP1_0018C" ~ "age_80_84",
           
           name == "DP1_0021C" ~ "age_18plus",
           name == "DP1_0024C" ~ "age_65_plus",
           
           name == "DP1_0073C" ~ "age_median",
           
           name == "DP1_0025C" ~ "male",
           name == "DP1_0049C" ~ "female",
           
           name == "DP1_0076C" ~ "race_total_pop",
           name == "DP1_0078C" ~ "race_white", #use?
           name == "DP1_0079C" ~ "race_black_aa",
           name == "DP1_0080C" ~ "race_native",
           name == "DP1_0081C" ~ "race_asian",
           name == "DP1_0082C" ~ "race_hawaiian_islander",
           name == "DP1_0083C" ~ "race_other",
           name == "DP1_0084C" ~ "race_multiple",
           
           name == "DP1_0093C" ~ "hispanic_total",
           name == "DP1_0094C" ~ "non_hispanic",
           name == "DP1_0105C" ~ "non_hispanic_white",
           
           name == "DP1_0115C" ~ "spouse_opposite_sex",
           name == "DP1_0116C" ~ "spouse_same_sex",
           name == "DP1_0117C" ~ "unmarried_opposite_sex",
           name == "DP1_0118C" ~ "unmarried_same_sex")) 


decennial_names_dp <- decennial_variables_dp %>% pull(name)

decennial_dp <- get_decennial(geography = "tract",
                              variables = decennial_names_dp,
                              year=yr,
                              sumfile="dp",
                              state = "WA") %>%
  filter(grepl(paste(study_counties, collapse = "|"), NAME)) %>%
  rename(variable0=variable) %>%
  left_join(select(decennial_variables_dp, name, variable), by=c("variable0"= "name"))  

########################################
# create new groups
new_demo_groups <- decennial_dp %>%
  select(GEOID, variable, value) %>%
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(
    # age groups
    age_00_19 = age_00_04+age_05_09+age_10_14+age_15_19,
    age_05_19 = age_05_09+age_10_14+age_15_19,
    age_20_64 = age_20_24+age_25_29+age_30_34+age_35_39+age_40_44+age_45_49+age_50_54+age_55_59+age_60_64,
    # race groups
    # race_total_pop is usually/always same as total_pop. may slightly over-esimate poc since nh white may be undercounted
    people_of_color = race_total_pop - non_hispanic_white,
    # partnered
    partner_same_sex = spouse_same_sex + unmarried_same_sex,
    partner_opposite_sex = spouse_opposite_sex + unmarried_opposite_sex
  ) %>%
  select(GEOID, age_00_19, age_05_19, age_20_64, people_of_color, partner_same_sex, partner_opposite_sex) %>%
  pivot_longer(cols = -GEOID, names_to = "variable", values_to = "value")

decennial_dp <- bind_rows(decennial_dp, new_demo_groups)

##################################################################################################
# housing 
decennial_variables_dhc <- load_variables(year = yr, dataset = "dhc") %>%
  filter(concept == "TENURE BY RACE OF HOUSEHOLDER",
         (label %in% c(" !!Total:", " !!Total:!!Owner occupied:", " !!Total:!!Renter occupied:"))) %>%
  mutate(dataset = "dhc",
         year = yr,
         variable = case_when(
           name == "H10_001N" ~ "house_total",
           name == "H10_002N" ~ "house_owners",
           name == "H10_010N" ~ "house_renters"))

decennial_names_dhc <- decennial_variables_dhc %>% pull(name)

decennial_dhc <- get_decennial(geography = "tract",
                               variables = decennial_names_dhc,
                               year=yr,
                               sumfile="dhc",
                               state = "WA") %>%
  filter(grepl(paste(study_counties, collapse = "|"), NAME)) %>%
  rename(variable0=variable) %>%
  left_join(select(decennial_variables_dhc, name, variable), by=c("variable0"= "name"))

##################################################################################################
# poverty

# ACS data differ from decennial Census data as they are based on an annual sample of approximately 3 million households, rather than a more complete enumeration of the US population. In turn, ACS data points are estimates characterized by a margin of error
# moe: margin of error (default 90% Conf lvl around estimate)

## every 5 yrs, block group data, 'estimates' with errors

# --> TO DO, add: education, language, citizenship/foreign born,ratio of income to poverty 
#--> see other options, e.g. "pl" https://walker-data.com/tidycensus/reference/load_variables.html

acs_variables_acs5 <- load_variables(year = yr, 
                                     dataset = "acs5" #cache=T
) %>%
  filter(geography == "tract",
         label %in% c("Estimate!!Total:!!Income in the past 12 months below poverty level:", "Estimate!!Total:!!Income in the past 12 months at or above poverty level:") & concept == "POVERTY STATUS IN THE PAST 12 MONTHS BY SEX BY AGE")  %>%
  mutate(dataset = "acs5",
         year = yr,
         variable = case_when(
           name == "B17001_002" ~ "poverty_below",
           name == "B17001_031" ~ "poverty_above")) %>%
  select(-geography)

acs_names_acs5 <- acs_variables_acs5 %>% pull(name)

# note that these are "estimates" not values
decennial_acs <- get_acs(geography = "tract",
                         variables = acs_names_acs5,
                         year=yr,
                         sumfile="acs5",
                         state = "WA") %>%
  filter(grepl(paste(study_counties, collapse = "|"), NAME)) %>%
  rename(variable0=variable) %>%
  left_join(select(acs_variables_acs5, name, variable), by=c("variable0"= "name")) %>%
  rename(value=estimate) %>%
  select(-c(moe))

##################################################################################################
# combine files 
pop_variables <- bind_rows(decennial_variables_dp,
                           decennial_variables_dhc,
                           acs_variables_acs5)

pop <- bind_rows(decennial_dp,
                 decennial_dhc,
                 decennial_acs) %>%
  filter(GEOID %in% tracts$GEOID) %>%
  select(-c(NAME)) %>%
  rename(pop=value)

##################################################################################################

#rm(api_key)

##################################################################################################
# VULNERABILITY INDICES
##################################################################################################
# CDC/ATSDR Social Vulnerability Index (SVI)
svi <- read.csv(file.path("data", "raw", "SVI_2020_US.csv")) %>%
  select(FIPS, LOCATION, contains("RPL_THEME")) %>%
  rename(GEOID = FIPS,
         RPL_THEME1_SES = RPL_THEME1,
         RPL_THEME2_HOUSING = RPL_THEME2,
         RPL_THEME3_RACE_ETH = RPL_THEME3,
         RPL_THEME4_HOUSING_TRANSPORT = RPL_THEME4,
         RPL_THEMES_OVERALL = RPL_THEMES) %>%
  mutate(GEOID = as.character(GEOID),
         index = "svi"
  ) %>%
  #filter(GEOID %in% tracts$GEOID) %>%
  left_join(select(tracts, GEOID), .) %>%
  st_drop_geometry() %>%
  select(-LOCATION) %>% 
  pivot_longer(cols=contains("THEME"), names_to = "index_variable", values_to = "index_value")

# higher values indicate more vulnerability
# ggplot() + geom_sf(data=svi, aes(fill=RPL_THEMES_OVERALL)) + labs(fill="Overall SVI")


##################################################################################################
# RISK ESTIMATES FROM STUDIES
##################################################################################################
# function returns the beta & se in log units
# note: function requires the RR and CI to be labeled: 'mean', 'lower', 'upper', and 'unit_increase'
get_beta_se_from_rr <- function(dt) {
  dt_names <- names(dt)
  
  # for a 95% CI, ~1.96
  z <- qnorm(.975) 
  
  # transform estimates to the native scale
  dt %>% 
    mutate(mean_log = log(mean), # beta
           # the log transformed confidence interval is symmetrical around the point estimate (could also use the delta method)
           ## log_confidence_bounds = log_mean +- 1.96*SE
           ## SE = (log_confidence_bouns - log_mean)/1.96  [SE = SD/sqrt(n)]
           ## SE is the average of these two intervals / 1.96
           lower_log = log(lower),
           upper_log = log(upper),
           se_log = ((upper_log - mean_log) + (mean_log - lower_log))/2/z,
           
           # transformed mean & SE to 1 unit
           mean_log_one_unit = mean_log/unit_increase,#*new_units,
           ## units of variability also change (e.g., https://stats.stackexchange.com/questions/241743/changing-units-of-measurement-in-a-simple-regression-model)
           se_log_one_unit = se_log/unit_increase) %>%
    select(all_of(dt_names), mean_log, se_log, mean_log_one_unit, se_log_one_unit)
}
##################################################################################################
risk <- read_xlsx(file.path("data", "raw", "ufp_health_outcomes.xlsx"), sheet = 1) %>%
  filter(include == "TRUE") %>%
  # estimate beta & SEs on the native scale and for a 1 unit change in exposure # check that all estimates are HR, RR, HR, etc.
  get_beta_se_from_rr()

##################################################################################################
# HEALTH IMPACT FUNCTIONS
##################################################################################################
# Sacks 2018: most users use a log-linear fn: deltaY = (1-EXP(-Beta*DELTAQ))*Incidence*POP
# see BenMap manual table 4-6 (or 4-9?)? see table 4-10 Beta Distribution Types and Variables
 
# WHO/EEA: "The main uncertainty is associated with the concentration-response functions used in the health risk assessment"
# use a CI to assess this uncertainty
# there are others cited regularly but for PM2.5 (e.g., Krewski https://www.healtheffects.org/system/files/Krewski140.pdf)
# RR is a 'log-linear' conc response fn: log(RR) = B*x

# --> TO DO: use alternatives? 

health_impact_functions <- data.frame(name = as.character(),
                                      form = as.character(),
                                      notes = as.character()) %>%
  add_row(name="attributable fraction approach", # with log-linear exposure response function
          form = "af*incidence*pop", 
          notes="European Environ Agency: https://www.eea.europa.eu/publications/assessing-the-risks-to-health;\nLiu 2020 wildfire paper")


##################################################################################################
# INCIDENCE RATES FOR EACH OUTCOME
##################################################################################################

# --> download more data from WA Tracking network? https://fortress.wa.gov/doh/wtn/WTNPortal/ 

# --> TO DO: denominator needs to be total births???

ptb <- data.frame(outcome = "PTB",
                  incidence = 0.093, 
                  population = "all",
                  geography = "county",
                  source = "King County: https://kingcounty.gov/en/legacy/depts/health/data/-/media/depts/health/data/documents/Health-of-Mothers-and-Infants-by-Race-Ethnicity.ashx#:~:text=Yearly%2C%20an%20average%20of%202%2C290,for%20a%20rate%20of%209.3%25."
)

##################################################################################################
# --> see EU health impact function that has been proposed
#--> could get incidence by subgroup (race/sex)?

# --> use total population as denominator? or just adults?

brain_tumor <- read.delim(file.path("data", "raw", "cdc_wonder", "United States and Puerto Rico Cancer Statistics, 1999-2020 Incidence.txt"), nrows = 1) %>%
  select(-(Notes)) %>%
  mutate(outcome = "brain tumor",
         incidence = Count/Population,
         population = "all",
         geography = "county",
         source = "CDC WONDER"
  ) %>%
  select(names(ptb))

##################################################################################################
# combine incidence data
incidence <- bind_rows(brain_tumor,
                       ptb)

##################################################################################################
# SAVE FILES
##################################################################################################
save(
  # aggregation units
  tracts,
  # population at aggregation units
  pop, pop_variables,
  # vulnerability indeces 
  svi,
  # air pollution at aggregation units
  air_pollution, 
  # health studies/risk 
  risk,
  health_impact_functions,
  # baseline incidence at aggregation units
  incidence,
  file = file.path(dt_path, "analysis_files.rda"))

