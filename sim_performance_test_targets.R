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
  tar_rep(sim_res_glmnet,
          command=run_targets_ltmle_simulation_comparison(
            library1="glmnet",
            SL.Control1=list(selector="undersmooth",alpha=1),
            library2="glmnet",
            SL.Control2=list(selector="undersmooth",alpha=0),
            n-n_df, 
            time=10,
            n_bootstrap_samples1=0,
            n_bootstrap_samples2=0,
            Markov_variables1=Markov_variables,
            Markov_variables2=Markov_variables),
          batches = 200, reps = 1, rep_workers = 1)

)


######################################################################
### realistic_targets.R ends here
