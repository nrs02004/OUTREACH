library("tidyr")
library("tidyverse")
library("lme4")
library("MASS")


## the formula we use to generate FEV (y) as a function of time (t) is y_{ik} = beta_{0i} + delta_i*t_{ik}+e_{ik}
## We actually have two models (given by the same formula, above, with different, but correlated parameters): One for clinic measured FEV and another for home measured FEV
## where observed variables are:
#### y_{ik} is the k-th obs for the i-th person
#### t_{ik} is the time of the k-th measurement on the i-th person (where t=1 is one year from baseline)
## parameters are:
#### beta_{0i} ~ N(beta_0_mean, beta_0_sd^2) : This is the mean FEV at time 0 for person i
#### delta_i ~ N(delta_mean, delta_sd^2) : This is the slope for FEV (from time 0 to 1) for person i
#### e_{ik} ~ N(0, residual_sd^2) : This is a noise variable (for each person at each time)
##### Again, the home vs clinic models have separate beta, delta and epsilon values --- the beta and delta values are correlated within person across those 2 models
## I believe you can think of delta_sd^2 as the between person variability for the change in FEV (from time 0 to 1)
## I believe you can think of residual_sd^2 as _related_to the within person variability for the change in FEV (it is not directly the variability)

### function gen_data_closure:
## Summary: takes in characteristics we would like to use for simulating the trial
### It then generates (and outputs) a function which can simulate a trial with a variable number of people

### Inputs:
## n_obs_per_person_clinic : The number of clinic observations for each person
## n_obs_per_person_home : The number of home observations for each person
## beta_0_mean_clinic : The overall mean for the person-specific random intercept in the clinic model
## beta_0_mean_home : The overall mean for the person-specific random intercept in the home model
## beta_0_sd_clinic : The standard deviation for the person-specific random intercept in the clinic model
## beta_0_sd_home : The standard deviation for the person-specific random intercept in the home model
## beta_0_cor : The correlation between beta_0 (the random intercept) for FEV measured at home and beta_0 for FEV measured in clinic
## delta_mean_clinic : The overall mean for the person-specific random slope in the clinic model
## delta_mean_home : The overall mean for the person-specific random slope in the home model
## delta_sd_clinic : The standard deviation for the person-specific random slope in the clinic model
## delta_sd_home : The standard deviation for the person-specific random slope in the home model
## delta_cor : The correlation between delta (the random slope) for FEV measured at home and delta for FEV measured in clinic
## residual_sd_clinic : The standard deviation for epsilon_{ik} in the clinic model
## residual_sd_home : The standard deviation for epsilon_{ik} in the home model
## followup_time : The length of the study in years

### Output:
## gen_data : a function with argument "n_people" that generates a trial with "n_people" using the model specified above

gen_data_closure <- function(n_obs_per_person_clinic,
                             n_obs_per_person_home,
                             beta_0_mean_clinic,
                             beta_0_mean_home,
                             beta_0_sd_clinic,
                             beta_0_sd_home,
                             beta_0_cor,
                             delta_mean_clinic,
                             delta_mean_home,
                             delta_sd_clinic,
                             delta_sd_home,
                             delta_cor,
                             residual_sd_clinic,
                             residual_sd_home,
                             followup_time){
  
  obs_times_clinic <- ((1:n_obs_per_person_clinic)-1)/(n_obs_per_person_clinic - 1)*followup_time
  obs_times_home <- ((1:n_obs_per_person_home)-1)/(n_obs_per_person_home - 1)*followup_time
  beta_0_cov = beta_0_sd_clinic * beta_0_sd_home * beta_0_cor
  delta_cov = delta_sd_clinic * delta_sd_home * delta_cor
  
  Sigma_beta_0 = matrix(c(beta_0_sd_clinic^2, beta_0_cov, beta_0_cov, beta_0_sd_home^2), ncol = 2)  
  Sigma_delta = matrix(c(delta_sd_clinic^2, delta_cov, delta_cov, delta_sd_home^2), ncol = 2)  
  
  beta_0_mean = c(beta_0_mean_clinic, beta_0_mean_home)
  delta_mean = c(delta_mean_clinic, delta_mean_home)
  
  gen_data <- function(n_people){

    beta_0 <- mvrnorm(n_people, mu = beta_0_mean, Sigma_beta_0)
    delta <- mvrnorm(n_people, mu = delta_mean, Sigma_delta)

    y_mean_clinic <-beta_0[,1] %o% rep(1,n_obs_per_person_clinic) + 
      delta[,1] %o% obs_times_clinic
    y_mean_home <-beta_0[,2] %o% rep(1,n_obs_per_person_home) + 
      delta[,2] %o% obs_times_home
    
    y_wide_clinic <- y_mean_clinic + rnorm(n_people*n_obs_per_person_clinic, 
                        sd = residual_sd_clinic)
    y_wide_clinic <- as.data.frame(y_wide_clinic)
    
    y_wide_home <- y_mean_home + rnorm(n_people*n_obs_per_person_home, 
                                           sd = residual_sd_home)
    y_wide_home <- as.data.frame(y_wide_home)
    
    colnames(y_wide_clinic) <- obs_times_clinic
    colnames(y_wide_home) <- obs_times_home
    y_wide_clinic[["ptid"]] <- rep(1:n_people)
    y_wide_home[["ptid"]] <- rep(1:n_people)
    y_long_clinic <- pivot_longer(y_wide_clinic, 
                           !ptid, 
                           names_to = "time",
                           values_to = "FEV") %>% mutate(loc = "clinic")
    y_long_home <- pivot_longer(y_wide_home, 
                                  !ptid, 
                                  names_to = "time",
                                  values_to = "FEV") %>% mutate(loc = "home")
    y_long = rbind(y_long_clinic, y_long_home)
    y_long$time = as.numeric(y_long$time)
    return(y_long)
  }
}

### function run_many:
## Summary: takes in characteristics of the data generating mechanism, and the analysis function...
## as well as the number of trials to simulate. It then summarizes the results of those simulations...
## in terms of vectors of widths corresponding to...
##    1) CI for bias of home vs clinic change in FEV; and
##    2) CI for the ratio of variances (or required sample sizes to obtain the same precision) for home vs clinic change in FEV 

### Inputs:
## data_gen_params : This is a list with all the input parameters to the "gen_data_closure" function
## ntrial : the number of trials to simulate
## nboot : the number of bootstrap replicates to use in our trial analysis
## alpha : the level corresponding to 1 - coverage of our two sided CI (used in our analysis)
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)

### Outputs:
## runs : a list containing output information from all the trials generated 
## mean_diff_width : the mean of the widths of each 1-alpha, 2-sided CI for the bias over the replicated datasets
## mean_diff_width : the mean of the widths of each 1-alpha, 2-sided CI for the ratio of variances over the replicated datasets


run_many <- function(data_gen_params,
                     n_people,
                     ntrial = 1e1,
                     nboot = 1e3,
                     alpha = 0.05,
                     only_2_clinic_measurements = FALSE){
  gen_dat <- gen_data_closure(data_gen_params$n_obs_per_person_clinic,
                              data_gen_params$n_obs_per_person_home,
                              data_gen_params$beta_0_mean_clinic,
                              data_gen_params$beta_0_mean_home,
                              data_gen_params$beta_0_sd_clinic,
                              data_gen_params$beta_0_sd_home,
                              data_gen_params$beta_0_cor,
                              data_gen_params$delta_mean_clinic,
                              data_gen_params$delta_mean_home,
                              data_gen_params$delta_sd_clinic,
                              data_gen_params$delta_sd_home,
                              data_gen_params$delta_cor,
                              data_gen_params$residual_sd_clinic,
                              data_gen_params$residual_sd_home,
                              data_gen_params$followup_time)
  runs <- replicate(ntrial, run_one(gen_dat, n_people, nboot, alpha, only_2_clinic_measurements))
  return(list(runs = runs,
              mean_diff_width = mean(unlist(runs[3,])),
              mean_ratio_width = mean(unlist(runs[4,]))))
}

### function run_one:
## Summary: This function is called by "run_many"; it generates data from a single trial...
## then analyzes that data and outputs a summary of that analysis

### Inputs:
## gen_dat : a function that generates trial data (and requires as input a number of observations)
## n_people : the number of people to include in the generated trial
## nboot : the number of bootstrap replicates to use in our trial analysis
## alpha : the level corresponding to 1 - coverage of our two sided CI (used in our analysis)
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)


### Outputs:
## out : a list containing
##       CI_diff : the confidence interval for the bias of change in FEV measured by home vs clinic spiro
##       CI_ratio : the confidence interval for the ratio of variances of estimated change in FEV measured by home vs clinic spiro
##       CI_diff_width : the width of CI_diff
##       CI_ratio_width : the width of CI_ratio
##       point_est : a 2-vector containing point estimates for the bias and ratio of variances

run_one <- function(gen_dat, n_people, nboot, alpha, only_2_clinic_measurements){
  dat <- gen_dat(n_people)
  dat_home <- dat %>% filter(loc == "home")
  dat_clinic <- dat %>% filter(loc == "clinic")
  out <- analyze_trial(dat_home, dat_clinic, nboot, alpha, only_2_clinic_measurements)
  return(out)
}

### function calc_stats:
## Summary: This function takes in data from a trial and outputs an estimate of...
## bias of change in FEV measured by home vs clinic spiro; and 
## the ratio of variances of change in FEV measured by home vs clinic spiro.

### Inputs:
## dat_home : a data frame containing home FEV measurements, times, and patient IDs
## dat_clinic : a data frame containing clinic FEV measurements, times, and patient IDs
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)

### Outputs:
## est_diff : estimated bias of change in FEV measured by home vs clinic spiro
## var_ratio : estimated ratio of variances of estimated change in FEV measured by home vs clinic spiro

calc_stats <- function(dat_home, dat_clinic, only_2_clinic_measurements){
  fit_home <- lmer(FEV~time + (1+time|ptid), data = dat_home)
  home_vals <- summary(fit_home)$coefficients[2,]
  if(only_2_clinic_measurements == FALSE){
    suppressMessages(
      fit_clinic <- lmer(FEV~time + (1+time|ptid), data = dat_clinic))
    clinic_vals <- summary(fit_clinic)$coefficients[2,]
  }else{
    suppressMessages(
      fit_clinic <- lm(FEV~-1+time + as.factor(ptid), data = dat_clinic))
    clinic_vals <- summary(fit_clinic)$coefficients[1,]
  }
  est_diff <- home_vals[1] - clinic_vals[1]
  var_ratio <- (home_vals[2] / clinic_vals[2])^2
  return(list(est_diff = est_diff, var_ratio = var_ratio))
}

### function analyze_trial:
## Summary: This function takes in data from a trial and uses a non-parametric bootstrap to output confidences for...
## the bias of change in FEV measured by home vs clinic spiro; and 
## the ratio of variances of estimated change in FEV measured by home vs clinic spiro.

### Inputs:
## dat_home : a data frame containing home FEV measurements, times, and patient IDs
## dat_clinic : a data frame containing clinic FEV measurements, times, and patient IDs
## nboot : the number of bootstrap replicates to use in calculating our CI
## alpha : the level corresponding to 1 - coverage of our two sided CI
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)


### Outputs:
## CI_diff : the confidence interval for the bias of change in FEV measured by home vs clinic spiro
## CI_ratio : the confidence interval for the ratio of variances of estimated change in FEV measured by home vs clinic spiro
## CI_diff_width : the width of CI_diff
## CI_ratio_width : the width of CI_ratio
## point_est : a 2-vector containing point estimates for the bias and ratio of variances

analyze_trial <- function(dat_home, dat_clinic, nboot, alpha, only_2_clinic_measurements){
  point_est <- calc_stats(dat_home, dat_clinic, only_2_clinic_measurements)
  boot_samps <- replicate(nboot, do_one_boot(dat_home, dat_clinic, only_2_clinic_measurements))
  
  CI_diff <- quantile(unlist(boot_samps[1,]), c(alpha/2, 1-alpha/2))
  CI_ratio <- quantile(unlist(boot_samps[2,]), c(alpha/2, 1-alpha/2))
  CI_diff_width = CI_diff[2] - CI_diff[1]
  CI_ratio_width = CI_ratio[2] - CI_ratio[1]
  return(list(CI_diff = CI_diff, CI_ratio = CI_ratio,
              CI_diff_width = CI_diff_width, CI_ratio_width = CI_ratio_width,
              point_est = point_est))
}
  
### function do_one_boot:
## Summary: Used by "analyze_trial" to create one bootstrapped replicated of our data...
## and calculate corresponding summaries on that bootstrapped dataset

### Inputs:
## dat_home : a data frame containing home FEV measurements, times, and patient IDs
## dat_clinic : a data frame containing clinic FEV measurements, times, and patient IDs
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)

### Outputs:
## stats : a list containing...
##         est_diff : estimated bias of change in FEV measured by home vs clinic spiro on bootstrapped data
##         var_ratio : estimated ratio of variances of estimated change in FEV measured by home vs clinic spiro on bootstrapped data


do_one_boot <- function(dat_home, dat_clinic, only_2_clinic_measurements){
  pts <- unique(dat_home$ptid)
  bootstrapped_id <- data.frame(ptid = sample(pts, replace = TRUE),
                                new_id = 1:length(pts))
  
  suppressMessages(resampled_dat_home <- left_join(bootstrapped_id, dat_home) %>% 
    dplyr::select(-ptid) %>% 
    rename(ptid = new_id))
  suppressMessages(resampled_dat_clinic <- left_join(bootstrapped_id, dat_clinic) %>% 
    dplyr::select(-ptid) %>% 
    rename(ptid = new_id))
  
  stats <- calc_stats(resampled_dat_home, resampled_dat_clinic, only_2_clinic_measurements)
  return(stats)
}

################ This is the end of the code used to generate/analyze simulated trials

################ The following functions are used to calculate the "true" ratio of variances (one of our parameters of interest) in our variance simulation scenarios

### function calc_ratio:
## Summary: This function uses simulation to identify/approximate...
## the true population ratio of variances of estimated change in FEV measured by home vs clinic spiro in a given simulation setting

### Inputs:
## data_gen_params : This is a list with all the input parameters to the "gen_data_closure" function
## n_people : the number of people to simulate in the trial: The ratio of variances should not depend on n (asymptotically, for sufficiently large n).
## ntimes : the number of trials to simulate
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)

### Outputs:
## ratio_calc : Our simulation-based approximation to the true population ratio of variances of estimated change in FEV measured by home vs clinic spiro

calc_ratio <- function(data_gen_params, n_people, ntimes = 1e3, only_2_clinic_measurements = TRUE){
  gen_dat <- gen_data_closure(data_gen_params$n_obs_per_person_clinic,
                              data_gen_params$n_obs_per_person_home,
                              data_gen_params$beta_0_mean_clinic,
                              data_gen_params$beta_0_mean_home,
                              data_gen_params$beta_0_sd_clinic,
                              data_gen_params$beta_0_sd_home,
                              data_gen_params$beta_0_cor,
                              data_gen_params$delta_mean_clinic,
                              data_gen_params$delta_mean_home,
                              data_gen_params$delta_sd_clinic,
                              data_gen_params$delta_sd_home,
                              data_gen_params$delta_cor,
                              data_gen_params$residual_sd_clinic,
                              data_gen_params$residual_sd_home)
  calcs <- replicate(ntimes, do_one_calc_ests(n_people, gen_dat, only_2_clinic_measurements))
  ests.home <- unlist(calcs[1,])
  ests.office <- unlist(calcs[2,])
  return(var(ests.home)/var(ests.office))
}

### function do_one_calc_ests:
## Summary: Helper function for "calc_ratio" --- this runs a single simulation, and calculates estimates of
## the variances of estimated change in FEV measured by home and clinic spiro.

### Inputs:
## n_people : the number of people to include in the generated trial
## gen_dat : a function that generates trial data (and requires as input a number of observations)
## only_2_clinic_measurements : a TRUE/FALSE flag indicating if we specified only 2 clinic measurements (in that case a fixed effects analysis is used)

### Outputs:
### out: a list containing
##      est_home : estimated change in FEV measured by home spiro
##      est_clinic : estimated change in FEV measured by clinic spiro

do_one_calc_ests <- function(n_people, gen_dat, only_2_clinic_measurements){
  dat <- gen_dat(n_people)
  dat_home <- dat %>% filter(loc == "home")
  dat_clinic <- dat %>% filter(loc == "clinic")
  fit_home <- lmer(FEV~time + (1+time|ptid), data = dat_home)
  home_vals <- summary(fit_home)$coefficients[2,]
  if(only_2_clinic_measurements == FALSE){
    suppressMessages(
      fit_clinic <- lmer(FEV~time + (1+time|ptid), data = dat_clinic))
    clinic_vals <- summary(fit_clinic)$coefficients[2,]
  }else{
    suppressMessages(
      fit_clinic <- lm(FEV~-1+time + as.factor(ptid), data = dat_clinic))
    clinic_vals <- summary(fit_clinic)$coefficients[1,]
  }
  return(list(est_home = home_vals[1], est_clinic = clinic_vals[1]))
}