
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


# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
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
seeds_null=sample(0:1000000, 500, replace=FALSE)



#system.time({res_glm_2= mclapply_targets_ltmle_simulation(null_sim=TRUE, time=2, n_cores=50, estimator="glm",seeds=seeds_null[1:10], library="glm")})

system.time({res_glm_markov_2= mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="glm markov",seeds=seeds_null, Markov_variables=Markov_variables, library="glm")})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))

system.time({res_glmnet_markov_2= mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="lasso markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_glmnet_undersmooth_markov_2= mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed lasso markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_ridge_markov_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="ridge markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_ridge_undersmooth_markov_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed ridge markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_EN_markov_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="EN markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_EN_undersmooth_markov_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed ridge markov",seeds=seeds_null, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))

# system.time({res_glm_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="glm",seeds=seeds_null, gbounds=c(0,1), library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# 
# system.time({res_glmnet_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="lasso",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# system.time({res_glmnet_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed lasso",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# system.time({res_ridge_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="ridge",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# load(paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# 
# cat("Hello word\n")
# system.time({res_ridge_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed ridge",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# cat("1: res_ridge_undersmooth_untruncated_2 run\n")
# 
# system.time({res_EN_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="EN",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
# system.time({res_EN_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed EN",seeds=seeds_null, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_glm_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="glm markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables, library="glm")})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_glmnet_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="lasso markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_glmnet_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed lasso markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_ridge_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="ridge markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_ridge_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed ridge markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_EN_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="EN markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))
system.time({res_EN_undersmooth_markov_untruncated_2=mclapply_targets_ltmle_simulation(null_sim=TRUE,n_cores=50, estimator="undersmoothed ridge markov",seeds=seeds_null, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})

save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seeds_null_v2.Rdata"))


