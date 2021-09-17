library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

sim.output <- read.csv("saved_files/simulation_outcomes.csv")
sim.params <- read.csv("saved_files/simulation_params.csv")

#load("saved_files/all_runs.Rdata")

#dim(runs_store[[1]][[1]][[1]][[1]])

head(sim.output)
head(sim.params)

sim.output %>% filter(n_home_obs == 12, cors == 0.0, n_ppl == 100)

## Calculates the halfwidth of a CI corresponding to a conf level of "power"...
## and adds it to the halfwidth of a 95% CI
## This gives us the cutoff (minus truth) with which we have ~ "power" power to reject at the 0.05 level
cutoff.calc <- function(width.of.95.CI, power){
  quant.use <- 1-(1-power)/2
  halfwidth <- qnorm(quant.use)/qnorm(0.95)*width.of.95.CI/2
  cutoff <- halfwidth + width.of.95.CI/2
  return(cutoff)
}

sim.output.update <- sim.output %>% 
  mutate(`pow0.8_bias_delta` = sapply(diff_mean_width, cutoff.calc, 0.8),
         `pow0.9_bias_delta` = sapply(diff_mean_width, cutoff.calc, 0.9),
         `pow0.95_bias_delta` = sapply(diff_mean_width, cutoff.calc, 0.95),
         `pow0.8_ratio_delta` = sapply(ratio_var_width, cutoff.calc, 0.8),
         `pow0.9_ratio_delta` = sapply(ratio_var_width, cutoff.calc, 0.9),
         `pow0.95_ratio_delta` = sapply(ratio_var_width, cutoff.calc, 0.95))

### I want the thresholds at which we have 80%/90% power, if the true ratio is 1. 
### We do this in a bit of a hackish way:
###### ideally we would rerun a simulation setting for 12 home obs that had var ratio ~ 1
###### In the absence of that I am going to interpolate the variance of our estimate between optimistic and neutral

ratio_use <- sim.output.update %>% filter(n_ppl == 100, cors == 0.4, n_home_obs == 12) %>% pull(true_ratio)
thresh.use <- sim.output.update %>% filter(n_ppl == 100, cors == 0.4, n_home_obs == 12) %>% dplyr::select(pow0.8_ratio_delta, pow0.9_ratio_delta, pow0.95_ratio_delta)
prop.through <- (1-ratio_use[2])/(ratio_use[3] - ratio_use[2])
thresh.interp <- (thresh.use[3,] - thresh.use[2,])*prop.through + thresh.use[2,]
thresh.interp + 1

### Doing the same thing with Bias (not strictly necessary, but it seems like we should probably be thinking about engaging with the same simulation scenario for both)

ratio_use <- sim.output.update %>% filter(n_ppl == 100, cors == 0.4, n_home_obs == 12) %>% pull(true_ratio)
thresh.use <- sim.output.update %>% filter(n_ppl == 100, cors == 0.4, n_home_obs == 12) %>% dplyr::select(pow0.8_bias_delta, pow0.9_bias_delta, pow0.95_bias_delta)
prop.through <- (1-ratio_use[2])/(ratio_use[3] - ratio_use[2])
thresh.interp <- (thresh.use[3,] - thresh.use[2,])*prop.through + thresh.use[2,]
thresh.interp
