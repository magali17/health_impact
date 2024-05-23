# Magali Blanco
# 05/02/2024
# script to run health impact asssessment 

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

pacman::p_load(tidyverse, sf, parallel)

in_path <- file.path("data", "modified")
out_path <- file.path("data", "output")

set.seed(1)
use_cores <-1

# run a subset of scenarios
test_run <- TRUE 

##################################################################################################
# DATA
##################################################################################################
load(file.path(in_path, "analysis_files.rda"))

##################################################################################################
# COMMON VARIABLES
##################################################################################################
main_health_impact_fn <- health_impact_functions %>% 
  filter(grepl("attributable fraction", name)) %>% 
  pull(form)

# percent reduction in air pollution (AP)
ap_controls <- c(1, seq(25, 100, 25))
# health outcomes we will evaluate 
outcomes <- unique(risk$outcome)
# census groups (e.g., age, sex, race...)
# groups of interest
# pop_groups <- unique(pop$variable)

## --> update
pop_groups <- c("pop_total", 
                "non_hispanic_white", "people_of_color",
                "male", "female",
                "age_00_04", "age_05_19", "age_20_64", "age_65_plus",
                "partner_same_sex", "partner_opposite_sex",
                "poverty_above", "poverty_below"
)



if(test_run == TRUE) {
  ap_controls <- c(25,75) 
  #outcomes <- outcomes[1]
  #pop_groups <- c("pop_total", "male", "female")
}
##################################################################################################
# FUNCTIONS 
##################################################################################################
# calcualte air pollution during control scenarios

control_scenarios <- function(dt=air_pollution,
                              control_rollback = "percent",
                              control_value) {
  dt %>%
    group_by(GEOID) %>%
    mutate(control_scenario = paste0(control_rollback, "_", control_value),
           control_value = ifelse(control_rollback == "percent", Values*(1-control_value/100),
                                  ifelse(control_rollback == "fixed", control_value, NA)),
           # note: Values must be greater than control_value, otherwise exp(b) (RR/OR/HR) is different
           delta_ap = Values - control_value,
           #delta_ap = control_value-Values 
    ) %>%
    ungroup()
}

##################################################################################################
# get RR/OR/HR for a specific unit change in air pollution

transform_rr <- function(dt, outcome.) {
  # for a 95% CI, ~1.96
  z <- qnorm(.975)
  
  # transform beta, SE, HR, and CI from a 1 unit change to a new/transformed unit change
  temp_rr <- filter(risk, outcome.==outcome)
  ## if multiple betas per outcome, take the average on the log-linear scale. OK?
  mean_log_one_unit <- mean(temp_rr$mean_log_one_unit)
  se_log_one_unit <- mean(temp_rr$se_log_one_unit)
  
  dt %>% 
    mutate(outcome = first(temp_rr$outcome), 
           mean_log_trans = mean_log_one_unit*delta_ap,
           ## units of variability also change (e.g., https://stats.stackexchange.com/questions/241743/changing-units-of-measurement-in-a-simple-regression-model)
           se_log_trans = se_log_one_unit*delta_ap, 
           # RR & CI for new units
           mean_trans = exp(mean_log_trans),
           lower_trans = exp(mean_log_trans - z*se_log_trans),
           upper_trans = exp(mean_log_trans + z*se_log_trans)) 
}

##################################################################################################
# calcualte attributable fraction for a given HR/OR/HR

attributable_fraction <- function(rr) {(rr-1)/rr}

##################################################################################################
# combine population/group of interest & outcome incidence information 
## (currently incidence info is for county-level estimates so does not change by group)

get_pop_info <- function(pop_variable, outcome.) {
  pop_temp <- filter(pop, variable %in% pop_variable)
  incidence_temp <- filter(incidence, outcome==outcome.)
  
  # --> NOTE: need to update code if have multiple incidence rows at some point
  if(incidence_temp$geography == "county" & incidence_temp$population == "all") {
    pop_temp <- pop_temp %>%
      mutate(outcome=outcome.,
             # average if multiple incidence estimates (unlikely?)
             incidence = mean(incidence_temp$incidence))
  }
  
  # [add other scenarios as needed]
  
  return(pop_temp)
}

##################################################################################################
# health impact function
## combines exposure, population, incidence information to estimate changes in an outcome related to a specific air pollution change for a specific group

calculate_a_health_impact <- function(exposure_dt, #has current & control scenario AP levels
                                      outcome_dt # has pop counts & area-level incidence info
                                      # --> TO DO: HOW TO PASS IN A STRING?
                                      #, health_function 
) {
  exposure_dt %>%
    select(GEOID, Values, control_scenario, control_value, delta_ap, outcome, contains("af_")) %>%
    pivot_longer(contains("af_"), names_to = "af_name", values_to = "af") %>%  
    left_join(outcome_dt, by=c("GEOID", "outcome")) %>% 
    
    # --> NEED TO UPDATE
    mutate(outcome_control = af*incidence*pop,
           #outcome_control = health_function
    )  
}

##################################################################################################
# function repeats calculate_a_health_impact() many times for multiple groups (e.g., age/sex/race), control scenarios, and outcomes of interest
health_impacts <- function(pop_groups.=pop_groups, outcomes.=outcomes, ap_controls.=ap_controls) {
  
  result <- mclapply(pop_groups., mc.cores = use_cores, function(this_pop_group) {
    # for a specific outcome
    lapply(outcomes., function(this_outcome){
      # for a specific control scenario
      lapply(ap_controls., function(this_control_scenario) {
        
        exposure_dt. <- air_pollution %>%
          # calculate exposure levels under control scenarios
          control_scenarios(dt=.,control_rollback = "percent", this_control_scenario) %>%
          # calculate a health beta for this AP exposure change
          transform_rr(dt = ., outcome. = this_outcome) %>% 
          # calculate attributable fraction based on a health HR
          mutate(af_mean = attributable_fraction(mean_trans),
                 af_lower = attributable_fraction(lower_trans),
                 af_upper = attributable_fraction(upper_trans))
        
        outcome_dt. <- get_pop_info(pop_variable = this_pop_group, outcome. = this_outcome)
        
        calculate_a_health_impact(exposure_dt = exposure_dt., outcome_dt = outcome_dt.)
      }) %>% bind_rows() 
    }) %>% bind_rows()
  })  %>% bind_rows()
  
  return(result)
}

##################################################################################################
# HEALTH IMPACT ASSESSMENT
##################################################################################################
# --> consider saving smaller grouped files? if(!file.exists() & override_file=TRUE)...

health_impact <- health_impacts(pop_groups.=pop_groups, outcomes.=outcomes, ap_controls.=ap_controls)

##################################################################################################
# SAVE RESULTS
##################################################################################################
save(health_impact, 
     pop_groups, #main variables to be used in (male/female, ages, race....)
     file = file.path(out_path, "health_impact_files.rda"))
