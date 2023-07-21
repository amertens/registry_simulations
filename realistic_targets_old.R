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

# Uncomment below to deploy targets to parallel jobs
# on a Sun Grid Engine cluster when running tar_make_clustermq().
# options(clustermq.scheduler = "sge", clustermq.template = "sge.tmpl")


library(crew)
tar_option_set(
  controller = crew_controller_local(workers = 64)
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


# # Set the number of simulations
 n_sims = 2
 
 #set time horizon
 time=2
 library="glm"
 
# 
# # Create a vector for the number of simulations
# sim_nums <- seq_len(n_sims)

# Create a function that contains your targets
simulation_targets <- function(){
   list(
    tar_target(simulated_data,
      simulate_data(lava_model = lava_model,n = 10000)
    )#,
    # tar_target(simulated_data_list,
    #            get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)
    # ),
    # tar_target(estimate,
    #            {v = run_ltmle(name_outcome="dementia",
    #                            time_horizon=time,
    #                            test=FALSE,
    #                            outcome_data=simulated_data_list$outcome_data,
    #                            regimen_data=simulated_data_list$regimen_data,
    #                            baseline_data=simulated_data_list$sim_baseline_covariates,
    #                            timevar_data=simulated_data_list$sim_time_covariates,
    #                            det.Q.function=NULL,# now build-in
    #                            SL.library=library,
    #                            #SL.cvControl=list(selector="undersmooth",alpha=1),
    #                            verbose=TRUE)
    #            }),
    # tar_target(table_estimate,{
    #   summary(estimate[[1]][[1]]$Ltmle_fit)
    # }
  )
  return(1)
}

#NOTE! Look into tar_rep_map to both pass parameters to the models and to do repetitions

# # Use tar_map to generate the targets across number of iterations
dynamic_targets <- tar_rep(
  name=sim_res,
  command=simulation_targets(),
  batches = 1,
  reps = 2,
  rep_workers = 2
)


targets_list <- list(
    tar_target(coefs,{
        source("data/coefs.txt")
        coefs
    }),
    # tar_target(
    #   index_batch,
    #   seq_len(2), # Change the number of simulation batches here.
    #   deployment = "main"
    # ),
    # tar_target(
    #   index_sim,
    #   seq_len(2), # Change the number of simulations per batch here.
    #   deployment = "main"
    # ),
    tar_target(lava_model,{
        get_lava_model(time_horizon = time, coefs = coefs)
    }),
    # tar_target(lava_model_A1,{
    #   get_lava_model(time_horizon = 10,coefs = coefs, counterfactual_A=1)
    # }),
    # tar_target(lava_model_A0,{
    #   get_lava_model(time_horizon = 10,coefs = coefs, counterfactual_A=0)
    # }),
    # tar_target(truth,{
    #   #update get_lava_model to calc truth
    #   #calc_truth(coefs, seed=12345, nsamp=1000000)
    # }),
    tar_target(sim_seeds, sample.int(n = 100000, size = 100, replace = FALSE))

)

target_list <- c(targets_list, dynamic_targets)
target_list

######################################################################
### realistic_targets.R ends here
