## Here I have some really ugly code used for generating output dataframes

widths.dt <- data.frame(scenarios = character(),
                        n_home_obs = numeric(),
                        cors = numeric(),
                        n_ppl = numeric(),
                        ratio_var_width = numeric(),
                        diff_mean_width = numeric(),
                        true_ratio = numeric(),
                        true_diff = numeric())
for(i in 1:length(scenarios)){
  for(j in 1:length(n_home_obs_per_person)){
    for(k in 1:length(cors)){
      for(l in 1:length(n_ppl_vec)){
        
        new_row = data.frame(scenarios = scenarios[[i]]$name,
                    n_home_obs = n_home_obs_per_person[j],
                    cors = cors[k],
                    n_ppl = n_ppl_vec[l],
                    diff_mean_width = widths[i,j,k,l,1],
                    ratio_var_width = widths[i,j,k,l,2],
                    true_ratio = true_ratio[i,j,k],
                    true_diff = 0)
        widths.dt <- rbind(widths.dt, new_row)
      }
    }
  }
}

true.params.dt <- data.frame(scenarios = character(),
                             n_home_obs = numeric(),
                             cors = numeric(),
                             true_ratio = numeric(),
                             true_diff = numeric())
for(i in 1:length(scenarios)){
  for(j in 1:length(n_home_obs_per_person)){
    for(k in 1:length(cors)){
      
      new_row = data.frame(scenarios = scenarios[[i]]$name,
                           n_home_obs = n_home_obs_per_person[j],
                           cors = cors[k],
                           true_ratio = true_ratio[i,j,k],
                           true_diff = 0)
      true.params.dt <- rbind(true.params.dt, new_row)
    }
  }
}

## The data frames stored here are in the form used for creating tables/plots in the rmd doc
write.csv(widths.dt, "saved_files/simulation_outcomes.csv", row.names = FALSE)
write.csv(true.params.dt, "saved_files/simulation_params.csv", row.names = FALSE)

### The info from all the runs is stored here (in case we need it for debugging) --- it is also used in generating sample CIs
save(list = c("runs_store"),file = "saved_files/all_runs.Rdata")

## Here we store things as R objects --- this was used previously to generate some figures (but I believe is no longer used)
save(list = c("widths","true_ratio"), file = "saved_files/output_for_tables_plots.Rdata")