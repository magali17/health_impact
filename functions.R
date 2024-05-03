
##################################################################################################
# ATTRIBUTABLE FRACTION
##################################################################################################
# attributable fraction
# --> correct? https://www.eea.europa.eu/publications/assessing-the-risks-to-health 

# attributable_fraction <- function(rr) {(rr-1)/rr}


##################################################################################################
# TRANSFORM RR/OR/HR FOR DIFFERENT UNITS
##################################################################################################
# estimate a HR & CI for a specific change in exposure (air pollution change)
# dt=air_pollution_control
# health_impact = health_impact_function
# exposure_change
# 
# transform_rr <- function(dt, health_impact, #exposure_change
#                          ) {
#   dt_names <- names(dt)
#   
#   # for a 95% CI, ~1.96
#   z <- qnorm(.975) 
#   
#   # # transform beta & SE from native scale & 1 unit to a HR with a new unit change
#    beta <- health_impact$mean_log
#    se_log <- health_impact$se_log
#     
#     dt %>% 
#       mutate(beta_trans = beta*delta_ap,
#              ## units of variability also change (e.g., https://stats.stackexchange.com/questions/241743/changing-units-of-measurement-in-a-simple-regression-model)
#              se_trans = se_log*delta_ap, 
#              # RR & CI for new units
#              mean_trans = exp(beta_trans),
#              lower_trans = exp(beta_trans - z*se_trans),
#              upper_trans = exp(beta_trans + z*se_trans)) %>%
#       select(all_of(dt_names), mean_trans, lower_trans, upper_trans)
#     }


##################################################################################################
# LOAD AERMOD FILES
##################################################################################################
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
#  
##################################################################################################




