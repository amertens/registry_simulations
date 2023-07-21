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

# Uncomment below to use local multicore computing
# when running tar_make_clustermq().
options(clustermq.scheduler = "multicore")


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
 n_sims = 2
 n_df=100000
 #set time horizon
 time=10
 SL.library="glm"
 Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                    "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                    "bb","ccb","rasi","thiazid",
                    "loop","mra","copd_med"  )
 
  
 #To do:
 #Make simple sim with just age and sex, no censoring or death
 #look into mclapply for parallelization (Thomas and Helena's project)
 #Set up bootstrap
 
 # hyperparameters <- tibble::tibble(
 #   #lava_model=list(lava_model,lava_model,lava_model),
 #   n=rep(10000,3),
 #   # library=c("glm","glmnet","glmnet"),
 #   # SL.cvControl=list(NULL,
 #   #                  list(selector="undersmooth",alpha=0),
 #   #                  list(selector="undersmooth",alpha=1))
 #   
 # )

list(
    # tar_target(coefs,{
    #     source("data/coefs.txt")
    #     coefs
    # }),
    # tar_target(lava_model,
    #     get_lava_model(time_horizon = time, coefs = coefs),
    #     deployment = "main"
    # ),
    # tar_target(res,{
    #   run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="undersmooth",alpha=1), n=n_df, time=time)
    # })
  
  tar_rep(sim_res_glm2,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=2),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm3,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=3),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm4,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=4),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm5,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=5),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm6,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=6),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm7,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=7),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm8,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=8),
          batches = 4, reps = 50, rep_workers = 50),
  tar_rep(sim_res_glm9,
          command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=9),
          batches = 4, reps = 50, rep_workers = 50),
    tar_rep(sim_res_glm,
            command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=time),
            batches = 4, reps = 50, rep_workers = 50),
    tar_rep(sim_res_glm_markov,
            command=run_targets_ltmle_simulation_batch(library="glm", n=n_df, time=time, Markov_variables=Markov_variables),
            batches = 4, reps = 50, rep_workers = 50),
    tar_rep(sim_res_glmnet,
            command=run_targets_ltmle_simulation_batch(library="glmnet", SL.Control=list(selector="min_lambda",alpha=1),  n=n_df, time=time),
            batches = 4, reps = 50, rep_workers = 50),
    tar_rep(sim_res_glmnet_undersmooth,
            command=run_targets_ltmle_simulation_batch(#library="glmnet", SL.Control=list(selector="undersmooth",alpha=1),
                                                       n=n_df, time=time),
            batches = 2, reps = 2, rep_workers = 2)
    ,tar_target(sim_res_tab_glm, clean_sim_res(res=sim_res_glm))
    ,tar_target(sim_res_tab_glmnet, clean_sim_res(res=sim_res_glmnet))
    ,tar_target(sim_res_tab_glmnet_undersmooth, clean_sim_res(res=sim_res_glmnet_undersmooth))
    #Set up to pass many hyperparameters to the simulation
    # tar_map_rep(sim_res,
    #         command=run_targets_ltmle_simulation(#lava_model=lava_model, 
    #                                              n=n, 
    #                                              library=library, time=time),
    #         values = hyperparameters,
    #           names = tidyselect::any_of("scenario"),
    #           batches = 2,
    #           reps = 3,
    #           rep_workers = 3
    #         )
)


######################################################################
### realistic_targets.R ends here
