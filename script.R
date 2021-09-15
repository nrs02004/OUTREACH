source("functions.R")
set.seed(1)

## Here we set vectors for some data generation parameters
### We loop through these vectors and simulate with each combination
n_ppl_vec <- c(50,75,100,125,150) ## Number of ppl in the study
n_home_obs_per_person = c(24, 12, 6) ## Number of home spiro measurements
cors <- c(0, 0.4, 0.8) ## Correlation of random effects between home and clinic spiro

followup_time <- 0.25 ## This is the length of the study in years; here we indicate the study will be 3 months

## We set some additional monte-carlo parameters
num_sims = 30 ## number of simulations per scenario/setting
nboot = 200 ## number of bootstrap replicates per simulation

## We now load in parameters for our 3 scenarios (optimistic, neutral, pessimistic)
### This creates the list "scenarios" with all our parameter values
source("scenarios.R")

## The following array will store average CI widths in all simulation settings; 
## widths[,,,,1] will contain the average CI widths for the bias parameter
## widths[,,,,2] will contain the average CI widths for the variance ratio parameter (the ratio of sample sizes that would be required with home vs clinic spiro to achieve the same power)

widths <- array(dim=c(length(scenarios), 
                      length(n_home_obs_per_person),
                      length(cors),
                      length(n_ppl_vec),
                      2))

## I use nested lists to store all the output from each of the runs...
## This seems much less elegent than storing things in a tibble where one column has entries that each contain all the runs from a specified simulation scenario
## Note. Check with Alex on how he does this!
runs_store <- list() 

## Now we loop through and run all the simulations
for(i in 1:length(scenarios)){
  runs_store[[i]] <- list()
  for(j in 1:length(n_home_obs_per_person)){
    runs_store[[i]][[j]] <- list()
    for(k in 1:length(cors)){
      runs_store[[i]][[j]][[k]] <- list()
      for(l in 1:length(n_ppl_vec)){
        params <- scenarios[[i]]
        params$beta_0_cor = cors[k]; params$delta_cor = cors[k]
        params$n_obs_per_person_home = n_home_obs_per_person[j]
        params$followup_time = followup_time
        runs <- run_many(data_gen_params = params,
                         n_people = n_ppl_vec[l],
                         ntrial = num_sims,
                         nboot = nboot,
                         alpha = 0.05,
                         only_2_clinic_measurements = (params$n_obs_per_person_home == 2))
        widths[i,j,k,l,1] <- runs$mean_diff_width
        widths[i,j,k,l,2] <- runs$mean_ratio_width
        runs_store[[i]][[j]][[k]][[l]] <- runs$runs
        print(paste(" ************* ", i,j,k,l ," *********** \n"))
      }
    }
  }
}

set.seed(2)
## Here we calculate/approximate the true population ratio of variances in each of our settings
## As a reminder, the true population bias is 0 in our simulation settings
true_ratio <- array(dim=c(length(scenarios), 
                          length(n_home_obs_per_person),
                          length(cors)))
for(i in 1:length(scenarios)){
  for(j in 1:length(n_home_obs_per_person)){
    for(k in 1:length(cors)){
      params <- scenarios[[i]]
      params$beta_0_cor = cors[k]; params$delta_cor = cors[k]
      params$n_obs_per_person_home = n_home_obs_per_person[j]
      params$followup_time = followup_time
      true_ratio[i,j,k] <- calc_ratio(data_gen_params=params, 
                                  n_people=150, 
                                  ntimes = 3e3,  
                                  only_2_clinic_measurements = FALSE)
    print(paste(" ************* ", i,j,k," *********** \n"))
    }
  }
}

true_ratio

## This script generates .csv and .Rdata files that will can be used later to summarize output
source("gen_output.R")
