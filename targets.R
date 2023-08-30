### realistic_targets.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 17 2023 (13:07) 
## Version: 
## Last-Updated: Jul 17 2023 (15:33) 
##           By: Thomas Alexander Gerds
##     Update #: 35
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
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

library(future)
library(future.callr)
plan(callr)

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
# Load LTMLE and augmentation functions
# avoid circular graph by adding the following line which sources all Ltmle functions
# to Ltmle.R and summary.Ltmle.R
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
 


list(
  #---------------------------------------------------------
  # Set up data simulation model and calculate truth
  #---------------------------------------------------------
    tar_target(coefs,{
        source("data/coefs.txt")
        coefs
    }),
    tar_target(lava_model,
        get_lava_model(time_horizon = time, coefs = coefs),
        deployment = "main"
    ),
  tar_target(seeds, sample(1:1000000, 500, replace=FALSE), deployment = "main"),
    tar_rep(truth_rep,
            command=calc_realistic_truth(),
            batches = 4, reps = 25, rep_workers = 25),
    tar_target(truth, average_truth(truth_rep)),

   tar_target(test, mclapply_targets_ltmle_simulation(estimator="test", seeds=seeds_1[1:2], library="glm",time=2))#,
  # tar_target(test2, mclapply_targets_ltmle_simulation(estimator="",seeds=seeds_1, library="glmnet",time=2, SL.Control=list(selector="undersmooth",alpha=1))),
  
  #---------------------------------------------------------
  # Run different estimation options
  #---------------------------------------------------------
  #  tar_target(res_glm_1, mclapply_targets_ltmle_simulation(estimator="glm",seeds=seeds_1, library="glm")),
  #  tar_target(res_glmnet_1, mclapply_targets_ltmle_simulation(estimator="lasso",seeds=seeds_1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))),
  #  tar_target(res_glmnet_undersmooth_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso",seeds=seeds_1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))),
  #  tar_target(res_ridge_1, mclapply_targets_ltmle_simulation(estimator="ridge",seeds=seeds_1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))),
  #  tar_target(res_ridge_undersmooth_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge",seeds=seeds_1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))),
  #  tar_target(res_EN_1, mclapply_targets_ltmle_simulation(estimator="EN",seeds=seeds_1, library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))),
  #  tar_target(res_EN_undersmooth_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed EN",seeds=seeds_1, library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))),
  # 
  # tar_target(res_glm_markov_1, mclapply_targets_ltmle_simulation(estimator="glm markov",seeds=seeds_1, Markov_variables=Markov_variables, library="glm")),
  # tar_target(res_glmnet_markov_1, mclapply_targets_ltmle_simulation(estimator="lasso markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))),
  # tar_target(res_glmnet_undersmooth_markov_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))),
  # tar_target(res_ridge_markov_1, mclapply_targets_ltmle_simulation(estimator="ridge markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))),
  # tar_target(res_ridge_undersmooth_markov_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))),
  # tar_target(res_EN_markov_1, mclapply_targets_ltmle_simulation(estimator="EN markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))),
  # tar_target(res_EN_undersmooth_markov_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds_1, Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))),
  # 
  # tar_target(res_glm_untruncated_1, mclapply_targets_ltmle_simulation(estimator="glm",seeds=seeds_1, gbounds=c(0,1), library="glm"))#,
  # tar_target(res_glmnet_untruncated_1, mclapply_targets_ltmle_simulation(estimator="lasso",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))),
  # tar_target(res_glmnet_undersmooth_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))),
  # tar_target(res_ridge_untruncated_1, mclapply_targets_ltmle_simulation(estimator="ridge",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))),
  # tar_target(res_ridge_undersmooth_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))),
  # tar_target(res_EN_untruncated_1, mclapply_targets_ltmle_simulation(estimator="EN",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))),
  # tar_target(res_EN_undersmooth_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed EN",seeds=seeds_1, gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))),
  # 
  # tar_target(res_glm_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="glm markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables, library="glm")),
  # tar_target(res_glmnet_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="lasso markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))),
  # tar_target(res_glmnet_undersmooth_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed lasso markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))),
  # tar_target(res_ridge_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="ridge markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))),
  # tar_target(res_ridge_undersmooth_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))),
  # tar_target(res_EN_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="EN markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))),
  # tar_target(res_EN_undersmooth_markov_untruncated_1, mclapply_targets_ltmle_simulation(estimator="undersmoothed ridge markov",seeds=seeds_1, gbounds=c(0,1), Markov_variables=Markov_variables,  library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))),

  #---------------------------------------------------------
  # Calculate simulation performance
  #---------------------------------------------------------
  # tar_target(sim_performance, calc_sim_performance(
  #   res=list(
  #     res_glm_1,
  #     res_glmnet_1,
  #     res_glmnet_undersmooth_1,
  #     res_ridge_1,
  #     res_ridge_undersmooth_1,
  #     res_EN_1,
  #     res_EN_undersmooth_1,
  #     res_glm_markov_1,
  #     res_glmnet_markov_1,
  #     res_glmnet_undersmooth_markov_1,
  #     res_ridge_markov_1,
  #     res_ridge_undersmooth_markov_1,
  #     res_EN_markov_1,
  #     res_EN_undersmooth_markov_1,
  #     res_glm_untruncated_1,
  #     res_glmnet_untruncated_1,
  #     res_glmnet_undersmooth_untruncated_1,
  #     res_ridge_untruncated_1,
  #     res_ridge_undersmooth_untruncated_1,
  #     res_EN_untruncated_1,
  #     res_EN_undersmooth_untruncated_1,
  #     res_glm_markov_untruncated_1,
  #     res_glmnet_markov_untruncated_1,
  #     res_glmnet_undersmooth_markov_untruncated_1,
  #     res_ridge_markov_untruncated_1,
  #     res_ridge_undersmooth_markov_untruncated_1,
  #     res_EN_markov_untruncated_1,
  #     res_EN_undersmooth_markov_untruncated_1
  #   ),
  #   truth=truth,
  #   time=10
  # ))
)


######################################################################
### realistic_targets.R ends here
