
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
seed5=sample(4000001:5000000, 200, replace=FALSE)

#run extra iteractions for estimation methods that have failed for some runs:
#                                 estimator N_reps
# 1     EN_undersmooth_markov_untruncated    468
# 2            EN_undersmooth_untruncated    181
# 3 glmnet_undersmooth_markov_untruncated    468
# 4        glmnet_undersmooth_untruncated    179
# 5  ridge_undersmooth_markov_untruncated    494
# 6         ridge_undersmooth_untruncated    430

load(paste0(here::here(),"/data/sim_results/sim_res_seed5.Rdata"))
res_EN_undersmooth_untruncated_5=res_EN_undersmooth_untruncated_4
rm(res_EN_undersmooth_untruncated_4)
# system.time({res_EN_undersmooth_untruncated_5=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed5, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed5.Rdata"))


system.time({res_glmnet_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed lasso",seeds=seed5[1:200], gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed5.Rdata"))

system.time({res_EN_undersmooth_untruncated_2=mclapply_targets_ltmle_simulation(n_cores=90, estimator="undersmoothed EN",seeds=seed5[1:200], gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_seed5.Rdata"))