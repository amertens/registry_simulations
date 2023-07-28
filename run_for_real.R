shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)

gc()
#tar_make(script = "realistic_targets.R")
tar_make(script = "null_targets.R")

sim_performance = tar_read(sim_performance)
sim_performance
write.csv(sim_performance, file=paste0(here::here(),"/data/sim_performance.csv"))

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


sim_res_glm_bootstrap = tar_read(sim_res_glm_bootstrap)
saveRDS(sim_res_glm_bootstrap, paste0(here::here(), "/scratch/sim_res_glm_bootstrap.RDS"))
