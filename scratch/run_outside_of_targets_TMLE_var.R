
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
set.seed(12345)
seeds=sample(0:1000000, 500, replace=FALSE)

res_glm_test_ic = mclapply_targets_ltmle_simulation(time=2, n_df=10000,
                                                estimator="glm",seeds=seeds[1:50], library="glm")
res_glm_test_tmle = mclapply_targets_ltmle_simulation(time=2, n_df=10000,
                                                    estimator="glm",seeds=seeds[1:50], library="glm",
                                                    tmle_var=TRUE)


#Run undersmoothed ridge
system.time({res_ridge_undersmooth_markov_tmle=mclapply_targets_ltmle_simulation(tmle_var=TRUE, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_tmle,  file=paste0(here::here(),"/data/sim_results/sim_res_undersmooth_ridge_markov_tmle.RDS"))

system.time({res_ridge_undersmooth_markov_null_tmle=mclapply_targets_ltmle_simulation(null_sim=TRUE, tmle_var=TRUE, n_cores=90, estimator="undersmoothed ridge markov",seeds=seeds, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
saveRDS(res_ridge_undersmooth_markov_null_tmle,  file=paste0(here::here(),"/data/sim_results/sim_res_undersmooth_ridge_markov_null_tmle.RDS"))


