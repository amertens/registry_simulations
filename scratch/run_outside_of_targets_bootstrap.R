
rm(list=ls())
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)
library(parallel)

gc()



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
set.seed(12345)
seeds_bootstrap<-sample(1:1000000, 500, replace=FALSE)

N_bootstraps=200
#test=mclapply_targets_ltmle_simulation(n_bootstrap_samples=2, time=2, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds1[1:2], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))

system.time({res_ridge_undersmooth_markov_boot1=mclapply_targets_ltmle_bootstrap_simulation(n_bootstrap_samples=N_bootstraps, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds_bootstrap[1:100], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_boot1, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot1.RDS"))

system.time({res_ridge_undersmooth_markov_boot2=mclapply_targets_ltmle_bootstrap_simulation(n_bootstrap_samples=N_bootstraps, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds_bootstrap[101:200], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_boot2, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot2.RDS"))

system.time({res_ridge_undersmooth_markov_boot3=mclapply_targets_ltmle_bootstrap_simulation(n_bootstrap_samples=N_bootstraps, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds_bootstrap[201:300], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_boot3, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot3.RDS"))

system.time({res_ridge_undersmooth_markov_boot4=mclapply_targets_ltmle_bootstrap_simulation(n_bootstrap_samples=N_bootstraps, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds_bootstrap[301:400], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_boot4, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot4.RDS"))

system.time({res_ridge_undersmooth_markov_boot5=mclapply_targets_ltmle_bootstrap_simulation(n_bootstrap_samples=N_bootstraps, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds_bootstrap[401:500], Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_boot5, file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot5.RDS"))

