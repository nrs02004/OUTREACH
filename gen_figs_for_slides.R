library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

sim.output <- read.csv("saved_files/simulation_outcomes.csv")
sim.params <- read.csv("saved_files/simulation_params.csv")

load("saved_files/all_runs.Rdata")

## Changing the name of the "neutral" scenario to "data driven"
sim.output$scenarios[which(sim.output$scenarios == "neutral")] <- "data driven"

### A plot for 24 home obs and cor of 0.4, evaluating how CI width of bias changes w/ sample size
### (can change number of home obs and cor --- this was a scenario of particular interest)
sim.output %>% filter(n_home_obs == 12, cors == 0.4) %>%
  ggplot(aes(x = n_ppl, y = diff_mean_width/2, color = scenarios)) +
  geom_line() + 
  xlab("Number of people enrolled") + ylab("95% CI half-width for Bias") + 
  scale_color_manual(values=c("#00ba38", "#619CFF", "#F8766D"))

ggsave("saved_files/bias_plot.pdf", width = 6, height = 4)

### A plot for 24 home obs and cor of 0.4, evaluating how CI width of variance-ratio changes w/ sample size
### (can change number of home obs and cor --- this was a scenario of particular interest)
sim.output %>% filter(n_home_obs == 12, cors == 0.4) %>% 
  ggplot(aes(x = n_ppl, y = ratio_var_width/2, color = scenarios)) +
  geom_line() + 
  xlab("Number of people enrolled") + ylab("95% CI half-width for Analytic Efficiency") + 
  scale_color_manual(values=c("#00ba38", "#619CFF", "#F8766D"))

ggsave("saved_files/SSratio_plot.pdf", width = 6, height = 4)


### Below, we do some pre-work here to generate pictures with representative CIs from various trials
## I'm sure there is a more compact way to do this (and it's terrible practice to have the numbers below hard coded)
## This was pretty done quick and dirty for design day (though it wasn't used). Probably not worth cleaning up, unless the figures will go in one of the formal documents


scenarios <- c("optimistic", "neutral", "pessimistic")
n_ppl_vec <- c(50,75,100,125,150) ## Number of ppl in the study
n_home_obs_per_person = c(24, 12, 6) ## Number of home spiro measurements
cors <- c(0, 0.4, 0.8)

CIs <- data.frame(scenario = character(),
                  n_home_obs_per_person = numeric(),
                  cor = numeric(),
                  n_ppl = numeric(),
                  type = character(),
                  lower = numeric(),
                  upper = numeric(),
                  trial = numeric())

for(i in 1:length(scenarios)){
  for(j in 1:length(n_home_obs_per_person)){
    for(k in 1:length(cors)){
      for(l in 1:length(n_ppl_vec)){
        CIs_bias <- runs_store[[i]][[j]][[k]][[l]][1,] %>%
          unlist() %>% 
          matrix(., ncol = 2, byrow = TRUE)
        CIs_ratio <- runs_store[[i]][[j]][[k]][[l]][2,] %>%
          unlist() %>% 
          matrix(., ncol = 2, byrow = TRUE)
        num = nrow(CIs_ratio)
        
        new_df = data.frame(lower = c(CIs_bias[,1], CIs_ratio[,1]))
        new_df$upper = c(CIs_bias[,2], CIs_ratio[,2])
        new_df$trial = c(1:num,1:num)
        new_df$type = c(rep("Bias", num), rep("Sample Size Ratio", num))
        
        
        
        new_df$scenario = scenarios[i]
        new_df$n_home_obs_per_person = n_home_obs_per_person[j]
        new_df$cor = cors[k]
        new_df$n_ppl = n_ppl_vec[l]
        
        CIs = rbind(CIs, new_df)
      }
    }
  }
}

colors = c(rgb(0/255,112/255,192/255,255/255), rgb(6/255,176/255,80/255,255/255), rgb(255/255,0,0,255/255))

### Here we specify the scenario to use for our plot (created below). These values can be changed to anything included in the simulat study
cor_use = 0.4
n_ppl_use = 100
scenario_use = "pessimistic"
n_home_obs_per_person_use = 12

scenario_ind = which(scenarios == scenario_use)

params.use <- sim.params %>% filter(n_home_obs == n_home_obs_per_person_use, cors == cor_use, scenarios == scenario_use)
params.use$true_ratio

## This creates a plot of representative CIs for bias (from our simulations) 
CIs %>% filter(cor == cor_use, n_ppl == n_ppl_use, scenario == scenario_use, n_home_obs_per_person == n_home_obs_per_person_use) %>%
  filter(type == "Bias") %>%
  ggplot(data = .,
         aes(x = trial, y = (upper+lower)/2, ymin = lower, ymax = upper)) + 
  geom_hline(data = tibble(type = c("Bias"),
                           yint = c(params.use$true_diff)),
             aes(yintercept = yint), linetype = "solid", size = 1,
             color = "#dd77dd") +
  geom_linerange(col = colors[scenario_ind]) +
  scale_x_continuous(name = "Simulated trial index") +
  ylim(-4,4) + 
  ylab("Bias CI")

ggsave(paste0("saved_files/bias_intervals_plot",
              "_scen=",scenario_use,
              "_n=",n_ppl_use,
              "_cor=",cor_use,
              "_obs-per=",n_home_obs_per_person_use,
              ".pdf"), width = 6, height = 4)

## This creates a plot of representative CIs for variance-ratio (from our simulations). Here I call it "sample size ratio" as the two are identical to first order (and we felt that was more intuitive)
CIs %>% filter(cor == cor_use, n_ppl == n_ppl_use, scenario == scenario_use, n_home_obs_per_person == n_home_obs_per_person_use) %>%
  filter(type == "Sample Size Ratio") %>%
  ggplot(data = .,
         aes(x = trial, y = (upper+lower)/2, ymin = lower, ymax = upper)) + 
  geom_hline(data = tibble(yint = c(params.use$true_ratio)),
             aes(yintercept = yint), linetype = "solid", size = 1,
             color = "#dd77dd") +
  geom_hline(data = tibble(yint = 1),
             aes(yintercept = yint), linetype = "solid", size = 1,
             color = rgb(0.7,0.7,0.7,0.7))+ 
  geom_linerange(col = colors[scenario_ind]) +
  scale_x_continuous(name = "Simulated trial index") +
  ylim(0,5) + 
  ylab("Sample Size Ratio CI")

ggsave(paste0("saved_files/ratio_intervals_plot",
              "_scen=",scenario_use,
              "_n=",n_ppl_use,
              "_cor=",cor_use,
              "_obs-per=",n_home_obs_per_person_use,
              ".pdf"), width = 6, height = 4)  

