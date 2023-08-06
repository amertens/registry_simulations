
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

cat("Running set 3:\n")


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
seed3=sample(2000001:3000000, 200, replace=FALSE)


load(paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glm_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm",seeds=seed3, library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# 
#  system.time({res_glmnet_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso",seeds=seed3, library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  system.time({res_glmnet_undersmooth_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso",seeds=seed3, library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  system.time({res_ridge_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge",seeds=seed3, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  system.time({res_ridge_undersmooth_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge",seeds=seed3, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  system.time({res_EN_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN",seeds=seed3, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  system.time({res_EN_undersmooth_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed3, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
#  save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
#  
# system.time({res_glm_markov_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm markov",seeds=seed3, Markov_variables=Markov_variables, library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glmnet_markov_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glmnet_undersmooth_markov_3= mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_ridge_markov_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_ridge_undersmooth_markov_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_EN_markov_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_EN_undersmooth_markov_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed3, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# 
# system.time({res_glm_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm",seeds=seed3, gbounds=c(0,1), library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glmnet_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
 system.time({res_glmnet_undersmooth_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
 save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_ridge_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_ridge_undersmooth_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_EN_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
 system.time({res_EN_undersmooth_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed3, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
 save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glm_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="glm markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables, library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glmnet_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="lasso markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_glmnet_undersmooth_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
system.time({res_ridge_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="ridge markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
 save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_ridge_undersmooth_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_EN_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="EN markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
# system.time({res_EN_undersmooth_markov_untruncated_3=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed ridge markov",seeds=seed3, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})

save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))