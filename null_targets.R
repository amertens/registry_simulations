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
tar_option_set(packages=c("fst","lava","ltmle","data.table","tidyverse","tictoc","glmnet","Matrix","Publish","matrixStats","speedglm","doParallel","parallel","caret","snow","doSNOW","foreach"))
tar_option_set(format = "qs")
tar_option_set(memory = "transient", garbage_collection = TRUE)
tar_option_set(storage = "worker", retrieval = "worker")
tar_option_set(error = "null")

# Uncomment below to use local multicore computing
# when running tar_make_clustermq().
options(clustermq.scheduler = "multicore")


#Set number of cores to parallelize the branches over
library(crew)
tar_option_set(
  controller = crew_controller_local(workers = 50)
)

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
# Load LTMLE and augmentation functions
# avoid circular graph by adding the following line which sources all Ltmle functions
# to Ltmle.R and summary.Ltmle.R
## nn=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)


# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
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

  #---------------------------------------------------------
  # Run different estimation options
  #---------------------------------------------------------
  #tar_target(test, run_targets_ltmle_simulation_null(library="glm", n=n_df, time=2)),
  tar_rep(null_res_glm,
          command=run_targets_ltmle_simulation_null(library="glm", n=n_df, time=time),
          batches = 4, reps = 50, rep_workers = 50),
          #batches = 1, reps = 1, rep_workers = 1),
  tar_rep(null_res_glm_markov,
          command=run_targets_ltmle_simulation_null(library="glm", n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(null_res_glmnet,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
          #batches = 1, reps = 1, rep_workers = 1)#,
  tar_rep(null_res_glmnet_undersmooth,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_glmnet_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_glmnet_undersmooth_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_EN,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_EN_undersmooth,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_EN_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_EN_undersmooth_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_ridge,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_ridge_undersmooth,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_ridge_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(null_res_ridge_undersmooth_markov,
          command=run_targets_ltmle_simulation_null(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  # tar_rep(null_res_RF,
  #         command=run_targets_ltmle_simulation_null(library="SL.randomForest",n=n_df, time=time),
  #         batches = 200, reps = 1, rep_workers = 1)
  
  tar_rep(sim_null_trunc_glm,
          command=run_targets_ltmle_simulation_trunc(library="glm", n=n_df, time=time),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_null_trunc_glm_markov,
          command=run_targets_ltmle_simulation_trunc(library="glm", n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_null_trunc_glmnet,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_glmnet_undersmooth,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_glmnet_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_glmnet_undersmooth_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_EN,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_EN_undersmooth,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_EN_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_EN_undersmooth_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_ridge,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_ridge_undersmooth,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_ridge_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_null_trunc_ridge_undersmooth_markov,
          command=run_targets_ltmle_simulation_trunc(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  #---------------------------------------------------------
  # Calculate results
  #---------------------------------------------------------

  tar_target(null_res_tab_glm, clean_sim_res(res=null_res_glm))
  ,tar_target(null_res_tab_glmnet, clean_sim_res(res=null_res_glmnet))
  ,tar_target(null_res_tab_glmnet_undersmooth, clean_sim_res(res=null_res_glmnet_undersmooth))
  ,tar_target(null_res_tab_glmnet_markov, clean_sim_res(res=null_res_glmnet_markov))
  ,tar_target(null_res_tab_glmnet_undersmooth_markov, clean_sim_res(res=null_res_glmnet_undersmooth_markov))
  ,tar_target(null_res_tab_ridge, clean_sim_res(res=null_res_ridge))
  ,tar_target(null_res_tab_ridge_undersmooth, clean_sim_res(res=null_res_ridge_undersmooth))
  ,tar_target(null_res_tab_ridge_markov, clean_sim_res(res=null_res_ridge_markov))
  ,tar_target(null_res_tab_ridge_undersmooth_markov, clean_sim_res(res=null_res_ridge_undersmooth_markov))
  ,tar_target(null_res_tab_EN, clean_sim_res(res=null_res_EN))
  ,tar_target(null_res_tab_EN_undersmooth, clean_sim_res(res=null_res_EN_undersmooth))
  ,tar_target(null_res_tab_EN_markov, clean_sim_res(res=null_res_EN_markov))
  ,tar_target(null_res_tab_EN_undersmooth_markov, clean_sim_res(res=null_res_EN_undersmooth_markov))
  #,tar_target(null_res_tab_RF, clean_sim_res(res=null_res_RF))
  
  #---------------------------------------------------------
  # Calculate simulation performance
  #---------------------------------------------------------
  
  ,tar_target(null_sim_performance, calc_sim_performance(
    res=list(
      null_res_tab_glm=null_res_tab_glm,
      null_res_tab_glmnet=null_res_tab_glmnet,
      null_res_tab_glmnet_undersmooth=null_res_tab_glmnet_undersmooth,
      null_res_tab_glmnet_markov=null_res_tab_glmnet_markov,
      null_res_tab_glmnet_undersmooth_markov=null_res_tab_glmnet_undersmooth_markov,
      null_res_tab_ridge=null_res_tab_ridge,
      null_res_tab_ridge_undersmooth=null_res_tab_ridge_undersmooth,
      null_res_tab_ridge_markov=null_res_tab_ridge_markov,
      null_res_tab_ridge_undersmooth_markov=null_res_tab_ridge_undersmooth_markov,
      null_res_tab_EN=null_res_tab_EN,
      null_res_tab_EN_undersmooth=null_res_tab_EN_undersmooth,
      null_res_tab_EN_markov=null_res_tab_EN_markov,
      null_res_tab_EN_undersmooth_markov=null_res_tab_EN_undersmooth_markov#,
      #null_res_tab_RF,null_res_tab_RF
    ),
    truth=data.frame(time=1:10, RR=rep(1,10), RD=rep(0,10)), 
    time=10
  ))
)


######################################################################
### realistic_targets.R ends here
