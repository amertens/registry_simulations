
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)
library(parallel)

gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)



# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)


# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
#set
#Longitudinal variables possibly following the markov process
Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                   "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                   "bb","ccb","rasi","thiazid",
                   "loop","mra","copd_med"  )

#seedlists
#seedlists
set.seed(12345)
seed1=sample(1:1000000, 500, replace=FALSE)

test=mclapply_targets_ltmle_simulation(n_bootstrap_samples=2, time=2, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds1[1:2], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))

#res_ridge_undersmooth_markov_boot=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))

#saveRDS(res_ridge_undersmooth_markov_boot, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot.RDS"))

