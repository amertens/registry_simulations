# JICI project
#
# Collaborators: Nerissa, Andrew, Bochra, Kathrine, Zeyi, Maya, Mark, Thomas, Christian
#
# Register data until 2018 are located at x:/Data/Rawdata_Hurtig/706582
# Corona update data are in the latest delivery (look for updates!)
# Currently we use
#  V:/Data/Workdata/706582/Corona_update/Data/11. levering.
#
# setwd("z:/Workdata/706582/Andrew Mertens/targets_diabetes_dementia/")
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"))
shh <-try(setwd("~/research/Methods/registry_simulations/"))
library(targets)
library(data.table)
library(tidyverse)
library(ltmle)
#load packages
tar_option_set(packages=c("lava","heaven","ltmle","data.table","tidyverse","SuperLearner","tictoc","glmnet","Matrix","Publish","matrixStats","speedglm","doParallel","parallel","caret","snow","doSNOW","foreach"))
tar_option_set(format = "qs")
tar_option_set(memory = "transient", garbage_collection = TRUE)
tar_option_set(storage = "worker", retrieval = "worker")

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nix1=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
#Load augmented LTMLE functions
 #nix2=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
nix3=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
nix3=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
 # source("./Ltmle/Augmentation/prepare_Ltmle.R")
 # source("./Ltmle/Augmentation/merge_data.R")
 # source("./Ltmle/Augmentation/get_subset_data.R")
 # source("./Ltmle/Augmentation/get_ltmle_data.R")
 # source("./Ltmle/Augmentation/get_formulas.R")
 # source("./Ltmle/Augmentation/get_rhs.R")

#To update:
  # #set up parallelization
  #  #use about half of space
  # ncores <- floor(detectCores()/2)
  # mycluster <- parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl=mycluster)


# -------------------------------------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------------------------

set.seed(12345)

#define some important global variables
dataset_N = 115698
sim_reps = 400

#baseline covariates
baseline_vars = c("ie_type","age_base","sex", "code5txt", "quartile_income")

#Longitudinal covariates
long_covariates = c("insulin_","any.malignancy_", "chronic.pulmonary.disease_","hypertension_",
                    "myocardial.infarction_", "ischemic.heart.disease_","heart.failure_", "renal.disease_", "sglt2_inhib_"   )

#Longitudinal covariates assumed to be following the markov process
Markov_variables = c("insulin","any.malignancy", "chronic.pulmonary.disease","hypertension",
                    "myocardial.infarction", "ischemic.heart.disease","heart.failure", "renal.disease", "sglt2_inhib")

# -------------------------------------------------------------------------------------------------------------
# Simulate data
# -------------------------------------------------------------------------------------------------------------

list(
  
#simulation parameters
tar_target(cc,fread(paste0(here::here(),"/data/coefficients.txt"))  #%>% 
  #drop censoring and death for simplicity
   # filter(!grepl("event_death",var) & !grepl("censor",var) ) 
)
,tar_target(sim_data, generate_data(cc, seed=12345, reps=sim_reps, n=dataset_N, N_time=10))

#calculate truth
,tar_target(truth,  calc_truth(cc, seed=12345, nsamp=1000000))




# -------------------------------------------------------------------------------------------------------------
# Run simulation
# -------------------------------------------------------------------------------------------------------------

#test single run
,tar_target(test_results_1, run_Ltmle(d=sim_data[[1]], SL.library = "glm", time_horizon=3, concurrentY=TRUE))
,tar_target(test_results_2, run_Ltmle(d=sim_data[[1]], SL.library = "glmnet", SL.cvControl=list(selector="min_lambda",alpha=1), time_horizon=3, concurrentY=TRUE))
,tar_target(test_results_3, run_Ltmle(d=sim_data[[1]], SL.library = "glmnet", SL.cvControl=list(selector="min_lambda",alpha=1), time_horizon=3, concurrentY=TRUE,  Markov_vars=Markov_variables))
#,tar_target(test_results_boot_2, run_ltmle_sim_bootstrap(sim_d_list=list(sim_data[[1]]), SL.library = "glm", time_horizon=3, concurrentY=TRUE, Niter=1, Nbootstrap=1))


# #test IC simulation with glm
# ,tar_target(glm_res, run_ltmle_sim(sim_d_list=sim_data, time_horizon=11, Ncores=64, Niter=200))
# #bootstrap test
# ,tar_target(test_results_bootstrap, run_ltmle_sim_bootstrap(sim_d_list=sim_data, 
#                                       SL.library = "glm",
#                                       #SL.cvControl=list(selector="min_lambda",alpha=1),
#                                       time_horizon=3, Ncores=64, Niter=2, Nbootstrap=2))
#IC simulation with glmnet
,tar_target(glm_res_t1, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glm",
                                         Markov_vars=NULL,
                                         time_horizon=2, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res_t1, run_ltmle_sim(sim_d_list=sim_data,
                                         SL.library = "glmnet",
                                         SL.cvControl=list(selector="min_lambda",alpha=1),
                                         Markov_vars=NULL,
                                         time_horizon=2, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res_t2, run_ltmle_sim(sim_d_list=sim_data,
                                         SL.library = "glmnet",
                                         SL.cvControl=list(selector="min_lambda",alpha=1),
                                         Markov_vars=NULL,
                                         time_horizon=3, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res_t4, run_ltmle_sim(sim_d_list=sim_data,
                                         SL.library = "glmnet",
                                         SL.cvControl=list(selector="min_lambda",alpha=1),
                                         time_horizon=5, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res_t6, run_ltmle_sim(sim_d_list=sim_data,
                                         SL.library = "glmnet",
                                         SL.cvControl=list(selector="min_lambda",alpha=1),
                                         time_horizon=7, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res_t8, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glmnet",
                                      SL.cvControl=list(selector="min_lambda",alpha=1),
                                      time_horizon=9, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glmnet",
                                      SL.cvControl=list(selector="min_lambda",alpha=1),
                                      time_horizon=11, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_undersmooth_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glmnet",
                                      SL.cvControl=list(selector="undersmooth",alpha=1),
                                      Markov_vars=NULL,
                                      time_horizon=11, Ncores=64, Niter=sim_reps))
,tar_target(glmnet_markov_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glmnet",
                                      SL.cvControl=list(selector="min_lambda",alpha=1),
                                      time_horizon=11, Ncores=64, Niter=sim_reps,
                                      Markov_vars=Markov_variables))
,tar_target(glmnet_markov_undersmooth_res, run_ltmle_sim(sim_d_list=sim_data,
                                                  SL.library = "glmnet",
                                                  SL.cvControl=list(selector="undersmooth",alpha=1),
                                                  time_horizon=11, Ncores=64, Niter=sim_reps,
                                                  Markov_vars=Markov_variables))
,tar_target(ridge_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glmnet",
                                      SL.cvControl=list(selector="min_lambda",alpha=0),
                                     Markov_vars=NULL,
                                      time_horizon=11, Ncores=64, Niter=sim_reps))

,tar_target(glm_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glm",
                                      time_horizon=11, Ncores=64, Niter=sim_reps))
,tar_target(glm_markov_res, run_ltmle_sim(sim_d_list=sim_data,
                                      SL.library = "glm",
                                      Markov_vars=Markov_variables,
                                      time_horizon=11, Ncores=64, Niter=sim_reps))

#bootstrap simulation with glmnet
# ,tar_target(glmnet_res_boot1, run_ltmle_sim_bootstrap(sim_d_list=sim_data,
#                                       SL.library = "glmnet",
#                                       SL.cvControl=list(selector="min_lambda",alpha=1),
#                                       time_horizon=11, Ncores=64, Niter=200, Nbootstrap=100))
# #IC simulation with glmnet and markov process
# ,tar_target(glmnet_res_markov, run_ltmle_sim(sim_d_list=sim_data,
#                                       SL.library = "glmnet",
#                                       SL.cvControl=list(selector="min_lambda",alpha=1),
#                                       time_horizon=11, Ncores=64, Niter=200, Markov_vars=Markov_variables))





#,tar_target(test_run,  run_ltmle(sim_data))

)
