### In here I have parameters for the simulations. 
## The variables stored in the lists below are explained in more depth at the top of "functions.R"
## They are passed into the function "gen_data_closure()" to create the data-generating mechanism

## The parameters in these models are for a year of followup (eg. delta_sd_clinic is the sd of the slope with units of years)

### This is a list of attributes that is shared among all simulation settings
shared = list(beta_0_mean_home = 80,
              beta_0_sd_home = 20,
              delta_mean_home = 0,
              n_obs_per_person_clinic = 3,
              beta_0_mean_clinic = 80,
              beta_0_sd_clinic = 20,
              delta_mean_clinic = 0,
              delta_sd_clinic = 2.6,
              residual_sd_clinic = 2.6)

### Below we define 3 specific settings (optimistic, neutral, and pessimistic)
## The variability of home measurements in each of these scenarios changes
optimistic <- list(delta_sd_home = 3, 
                   residual_sd_home = 3,
                   name = "optimistic")

neutral  <- list(delta_sd_home = 3.7, 
                 residual_sd_home = 3.7,
                 name = "neutral")

pessimistic  <- list(delta_sd_home = 5.2, 
                     residual_sd_home = 5.2,
                     name = "pessimistic")

### the correlation between home and clinic spiro measurements is not set here...
## A number of values are considered, and those are defined in "script.R"

scenarios <- list()
scenarios[[1]] <- optimistic
scenarios[[2]] <- neutral
scenarios[[3]] <- pessimistic

for(i in 1:length(scenarios)){
  scenarios[[i]] <- c(scenarios[[i]],shared)
}

