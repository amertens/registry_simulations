shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)

gc()
tar_make(script = "realistic_targets.R")

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
