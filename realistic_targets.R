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
    tar_rep(truth_rep,
            command=calc_realistic_truth(),
            batches = 4, reps = 25, rep_workers = 25),
    tar_target(truth, average_truth(truth_rep)),
    # tar_target(test,{
    #   run_targets_ltmle_simulation_bootstrap(library="glm",  n=n_df, time=2, n_bootstrap_samples=2)
    # }),
    # tar_rep(test2,
    #         command=run_targets_ltmle_simulation_bootstrap(library="glm",  n=n_df, time=2, n_bootstrap_samples=2),
    #         batches = 2, reps = 1, rep_workers = 1),
    # tar_target(test_tab, clean_sim_res(res=test)),

  
  #---------------------------------------------------------
  # Run different estimation options
  #---------------------------------------------------------
  tar_rep(sim_res_glm2,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=2),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm3,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=3),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm4,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=4),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm5,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=5),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm6,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=6),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm7,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=7),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm8,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=8),
          batches = 4, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm9,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=9),
          batches = 4, reps = 50, rep_workers = 25),
    tar_rep(sim_res_glm,
            command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=time),
            batches = 4, reps = 50, rep_workers = 25),
    tar_rep(sim_res_glm_markov,
            command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=time, Markov_variables=Markov_variables),
            batches = 4, reps = 50, rep_workers = 25),
    tar_rep(sim_res_glmnet,
            command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time),
            batches = 200, reps = 1, rep_workers = 1),
    tar_rep(sim_res_glmnet_undersmooth,
            command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time),
            batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_glmnet_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_glmnet_undersmooth_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_EN,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_EN_undersmooth,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_EN_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_EN_undersmooth_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_ridge,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_ridge_undersmooth,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_ridge_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=0),  n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_ridge_undersmooth_markov,
          command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=0),n=n_df, time=time, Markov_variables=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1),
  # tar_rep(sim_res_RF,
  #         command=run_targets_ltmle_simulation_batch(library="SL.randomForest",n=n_df, time=time),
  #         batches = 200, reps = 1, rep_workers = 1),
  tar_rep(sim_res_glm_1000iter,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=time),
          batches = 20, reps = 50, rep_workers = 25),
  tar_rep(sim_res_glm_bootstrap,
          command=run_targets_ltmle_simulation_bootstrap(library="glm",  n=n_df, time=10, n_bootstrap_samples=200),
          batches = 200, reps = 1, rep_workers = 1),
  #---------------------------------------------------------
  # Calculate results
  #---------------------------------------------------------
   tar_target(sim_res_tab_glm2, clean_sim_res(res=sim_res_glm2))
  ,tar_target(sim_res_tab_glm3, clean_sim_res(res=sim_res_glm3))
  ,tar_target(sim_res_tab_glm4, clean_sim_res(res=sim_res_glm4))
  ,tar_target(sim_res_tab_glm5, clean_sim_res(res=sim_res_glm5))
  ,tar_target(sim_res_tab_glm6, clean_sim_res(res=sim_res_glm6))
  ,tar_target(sim_res_tab_glm7, clean_sim_res(res=sim_res_glm7))
  ,tar_target(sim_res_tab_glm8, clean_sim_res(res=sim_res_glm8))
  ,tar_target(sim_res_tab_glm9, clean_sim_res(res=sim_res_glm9))
  ,tar_target(sim_res_tab_glm, clean_sim_res(res=sim_res_glm))
  ,tar_target(sim_res_tab_glmnet, clean_sim_res(res=sim_res_glmnet))
  ,tar_target(sim_res_tab_glmnet_undersmooth, clean_sim_res(res=sim_res_glmnet_undersmooth))
  ,tar_target(sim_res_tab_glmnet_markov, clean_sim_res(res=sim_res_glmnet_markov))
  ,tar_target(sim_res_tab_glmnet_undersmooth_markov, clean_sim_res(res=sim_res_glmnet_undersmooth_markov))
  ,tar_target(sim_res_tab_ridge, clean_sim_res(res=sim_res_ridge))
  ,tar_target(sim_res_tab_ridge_undersmooth, clean_sim_res(res=sim_res_ridge_undersmooth))
  ,tar_target(sim_res_tab_ridge_markov, clean_sim_res(res=sim_res_ridge_markov))
  ,tar_target(sim_res_tab_ridge_undersmooth_markov, clean_sim_res(res=sim_res_ridge_undersmooth_markov))
  ,tar_target(sim_res_tab_EN, clean_sim_res(res=sim_res_EN))
  ,tar_target(sim_res_tab_EN_undersmooth, clean_sim_res(res=sim_res_EN_undersmooth))
  ,tar_target(sim_res_tab_EN_markov, clean_sim_res(res=sim_res_EN_markov))
  ,tar_target(sim_res_tab_EN_undersmooth_markov, clean_sim_res(res=sim_res_EN_undersmooth_markov))
  #,tar_target(sim_res_tab_RF, clean_sim_res(res=sim_res_RF))
  
  #---------------------------------------------------------
  # Calculate simulation performance
  #---------------------------------------------------------
  
  ,tar_target(sim_performance, calc_sim_performance(
    res=list(
      sim_res_tab_glm=sim_res_tab_glm,
      sim_res_tab_glmnet=sim_res_tab_glmnet,
      sim_res_tab_glmnet_undersmooth=sim_res_tab_glmnet_undersmooth,
      sim_res_tab_glmnet_markov=sim_res_tab_glmnet_markov,
      sim_res_tab_glmnet_undersmooth_markov=sim_res_tab_glmnet_undersmooth_markov,
      sim_res_tab_ridge=sim_res_tab_ridge,
      sim_res_tab_ridge_undersmooth=sim_res_tab_ridge_undersmooth,
      sim_res_tab_ridge_markov=sim_res_tab_ridge_markov,
      sim_res_tab_ridge_undersmooth_markov=sim_res_tab_ridge_undersmooth_markov,
      sim_res_tab_EN=sim_res_tab_EN,
      sim_res_tab_EN_undersmooth=sim_res_tab_EN_undersmooth,
      sim_res_tab_EN_markov=sim_res_tab_EN_markov,
      sim_res_tab_EN_undersmooth_markov=sim_res_tab_EN_undersmooth_markov#,
      #sim_res_tab_RF,sim_res_tab_RF
    ), 
    truth=truth, 
    time=10
  ))
)


######################################################################
### realistic_targets.R ends here
