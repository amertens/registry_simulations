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
try(setwd("~/research/Methods/registry_simulations/"))
library(targets)
library(data.table)
library(tidyverse)
#load packages
tar_option_set(packages=c("lava","heaven","ltmle","data.table","tidyverse","SuperLearner","tictoc","glmnet","Matrix","Publish","matrixStats","speedglm","doParallel","parallel","caret","snow","doSNOW","foreach")
               )

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nix=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
source("Ltmle/Augmentation/Ltmle.R")
  # #set up parallelization
  #  #use about half of space
  # ncores <- floor(detectCores()/2)
  # mycluster <- parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl=mycluster)

# -------------------------------------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------------------------

set.seed(4534)

#define some important global variables
dataset_N = 1000
sim_reps = 10

#baseline covariates
baseline_vars = c("ie_type","age_base","sex", "code5txt", "quartile_income")

#Longitudinal covariates
long_covariates = c("insulin_","any.malignancy_", "chronic.pulmonary.disease_","hypertension_",
                    "myocardial.infarction_", "ischemic.heart.disease_","heart.failure_", "renal.disease_", "sglt2_inhib_"   )


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
,tar_target(truth,  calc_truth(cc, seed=23426, nsamp=100000))
)



# -------------------------------------------------------------------------------------------------------------
# Run simulation
# -------------------------------------------------------------------------------------------------------------


