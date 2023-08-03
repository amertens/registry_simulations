
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
tar_option_set(packages=c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","Publish","matrixStats","speedglm","parallel","caret","foreach","clustermq"))
tar_option_set(format = "qs")
tar_option_set(memory = "transient", garbage_collection = TRUE)
tar_option_set(error = "null")


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
seeds1=tar_read(seeds_1)


res_glm_1= mclapply_targets_ltmle_simulation(estimator="glm",seeds=seeds1, library="glm")
system.time({res_glmnet_1=mclapply_targets_ltmle_simulation(estimator="lasso",seeds=seeds1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})


#  res_glmnet_undersmooth_1= mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso",seeds=seeds1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))
#  res_ridge_1= mclapply_targets_ltmle_simulation(estimator="ridge",seeds=seeds1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))
#  res_ridge_undersmooth_1= mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge",seeds=seeds1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))
#  res_EN_1= mclapply_targets_ltmle_simulation(estimator="EN",seeds=seeds1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))
#  res_EN_undersmooth_1= mclapply_targets_ltmle_simulation(estimator="undersmoothed EN",seeds=seeds1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))
# 
# res_glm_markov_1= mclapply_targets_ltmle_simulation(estimator="glm markov",seeds=seeds1, Markov_variables=Markov_variables, library="glm")
# res_glmnet_markov_1= mclapply_targets_ltmle_simulation(estimator="lasso markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))
# res_glmnet_undersmooth_markov_1= mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))
# res_ridge_markov_1=mclapply_targets_ltmle_simulation(estimator="ridge markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))
# res_ridge_undersmooth_markov_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))
# res_EN_markov_1=mclapply_targets_ltmle_simulation(estimator="EN markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))
# res_EN_undersmooth_markov_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))
# 
# res_glm_untruncated_1=mclapply_targets_ltmle_simulation(estimator="glm",seeds=seeds1, gbounds=c(0,1), library="glm")
# res_glmnet_untruncated_1=mclapply_targets_ltmle_simulation(estimator="lasso",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))


res_glmnet_undersmooth_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))
res_ridge_untruncated_1=mclapply_targets_ltmle_simulation(estimator="ridge",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))
res_ridge_undersmooth_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))
res_EN_untruncated_1=mclapply_targets_ltmle_simulation(estimator="EN",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))
res_EN_undersmooth_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed EN",seeds=seeds1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))

res_glm_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="glm markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables, library="glm")
res_glmnet_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="lasso markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))
res_glmnet_undersmooth_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))
res_ridge_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="ridge markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))
res_ridge_undersmooth_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))
res_EN_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="EN markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))
res_EN_undersmooth_markov_untruncated_1=mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))


save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds1.Rdata"))

