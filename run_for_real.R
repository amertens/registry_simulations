shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)

gc()
tar_make(script = "realistic_targets.R")


test=run_targets_ltmle_simulation_comparison(rep_seed=12345, library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=10000, time=2, Markov_variables=Markov_variables)

#tar_make(script = "null_targets.R")
# tar_make(script = "sim_performance_test_targets.R")
# 
# n_df=100000
# test=run_targets_ltmle_simulation_comparison(
#   library1="glmnet",
#   SL.Control1=list(selector="undersmooth",alpha=1),
#   library2="glmnet",
#   SL.Control2=list(selector="undersmooth",alpha=0),
#   n-n_df, 
#   time=10,
#   n_bootstrap_samples1=0,
#   n_bootstrap_samples2=0,
#   Markov_variables1=Markov_variables,
#   Markov_variables2=Markov_variables)

# sim_performance = tar_read(sim_performance)
# sim_performance
# write.csv(sim_performance, file=paste0(here::here(),"/data/sim_performance.csv"))

# null_sim_performance  = tar_read(null_sim_performance )
# null_sim_performance 
# write.csv(null_sim_performance , file=paste0(here::here(),"/data/null_sim_performance .csv"))

#Need to debug, doesn't load functions
#tar_make_clustermq()

# #See how long each target ran:
# tar_visnetwork(script = "realistic_targets.R", label = c("time", "size", "branches"))
# 
# 
# # Shiny app to monitor progress
# tar_watch(seconds = 10, outdated = FALSE, targets_only = TRUE)
# # Now run the pipeline and watch the graph change.
# px <- tar_make()


# sim_res_glm_bootstrap = tar_read(sim_res_glm_bootstrap)
# saveRDS(sim_res_glm_bootstrap, paste0(here::here(), "/scratch/sim_res_glm_bootstrap.RDS"))
