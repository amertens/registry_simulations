shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(targets)
tar_make(script = "realistic_targets.R")

tar_make(simulated_data_list,script = "realistic_targets.R")
