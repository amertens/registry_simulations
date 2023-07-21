shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)
tar_make(script = "realistic_targets.R")

#Need to debug, doesn't load functions
#tar_make_clustermq()

#See how long each target ran:
tar_visnetwork(script = "realistic_targets.R", label = c("time", "size", "branches"))


# Shiny app to monitor
tar_watch(seconds = 10, outdated = FALSE, targets_only = TRUE)
# Now run the pipeline and watch the graph change.
px <- tar_make()

#TEMPORARY: check truth. Need to make function that averages across many seeds
counterfactual_data_A0  <- clean_sim_data(tar_read(counterfactual_data_A0), N_time=10)
counterfactual_data_A1  <- clean_sim_data(tar_read(counterfactual_data_A1), N_time=10)


table(counterfactual_data_A0$GLP1RA_9)
table(counterfactual_data_A1$GLP1RA_9)

mean(counterfactual_data_A0$dementia_10)*100
mean(counterfactual_data_A1$dementia_10)*100
mean(counterfactual_data_A0$dementia_10)*100-mean(counterfactual_data_A1$dementia_10)*100
mean(counterfactual_data_A0$dementia_10)-mean(counterfactual_data_A1$dementia_10)


truth <- tar_read(truth)
res <- tar_read(sim_res_glm2)
res_tab=NULL
for(i in 1:(length(res)/4)){
  temp = summary(res[[i + (3*(i-1))]]$GLP1RA$Ltmle_fit)
  res_tab <- bind_rows(res_tab, temp)
}