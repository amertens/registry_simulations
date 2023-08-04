
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})



#Set number of cores to parallelize the branches over

# library(crew)
# tar_option_set(
#   controller = crew_controller_local(workers = 50)
# )

#debugging parallelization issue
#https://github.com/ropensci/drake/issues/675#issuecomment-458222414

# Uncomment below to use local multicore computing
# when running tar_make_clustermq().
# tar_option_set(storage = "worker", retrieval = "worker", error="continue")
# options(clustermq.scheduler="multiprocess")



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
seed2=sample(1000001:2000000, 250, replace=FALSE)


# system.time({res_glm_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm",seeds=seed2, library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# 
#  system.time({res_glmnet_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso",seeds=seed2, library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  system.time({res_glmnet_undersmooth_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso",seeds=seed2, library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  system.time({res_ridge_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge",seeds=seed2, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  system.time({res_ridge_undersmooth_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge",seeds=seed2, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  system.time({res_EN_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN",seeds=seed2, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  system.time({res_EN_undersmooth_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed2, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
#  
# system.time({res_glm_markov_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm markov",seeds=seed2, Markov_variables=Markov_variables, library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_glmnet_markov_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_glmnet_undersmooth_markov_2= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_ridge_markov_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_ridge_undersmooth_markov_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_EN_markov_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_EN_undersmooth_markov_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed2, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# 
# system.time({res_glm_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm",seeds=seed2, gbounds=c(0,1), library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_glmnet_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_glmnet_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
# system.time({res_ridge_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))


system.time({res_ridge_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_EN_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_EN_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed2, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_glm_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables, library="glm")})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_glmnet_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_glmnet_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_ridge_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_ridge_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_EN_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
system.time({res_EN_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed2, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})

save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))

source(paste0(here::here(),"/scratch/run_outside_of_targets_3.R"))

