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
tar_option_set(packages=c("lava","ltmle","data.table","tidyverse","tictoc","glmnet","Matrix","Publish","matrixStats","speedglm","doParallel","parallel","caret","snow","doSNOW","foreach"))
tar_option_set(format = "qs")
tar_option_set(memory = "transient", garbage_collection = TRUE)
tar_option_set(storage = "worker", retrieval = "worker")
# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
# Load LTMLE and augmentation functions
# avoid circular graph by adding the following line which sources all Ltmle functions
# to Ltmle.R and summary.Ltmle.R
## nn=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

list(
    tar_target(coefs,{
        source("data/coefs.txt")
        coefs
    }),
    tar_target(lava_model,{
        get_lava_model(time_horizon = 10,coefs = coefs)
    }),
    tar_target(simulated_data,{
        simulate_data(lava_model = lava_model,n = 100000)
    }),
    tar_target(simulated_data_list,{
        get_simulated_data_list(simulated_data = simulated_data)
    }),
    tar_target(elastic_estimate,
    {
        ## tar_load(simulated_data_list)
        ## PL=prepare_Ltmle(regimen_data=simulated_data_list$regimen_data[[1]],outcome_data=simulated_data_list$outcome_data,name_outcome="dementia",name_regimen="GLP1RA",name_censoring = "Censored",censored_label = "0",name_comp.event = "Dead",baseline_data=simulated_data_list$sim_baseline_covariates,timevar_data=simulated_data_list$sim_time_covariates[,1:11,with = FALSE],time_horizon=10,SL.library="glm",Markov="",deterministic.Q.function=NULL,abar=list(rep(1,10),rep(0,10)))
        v = run_ltmle(name_outcome="dementia",
                      time_horizon=c(10),
                      test=FALSE,
                      outcome_data=simulated_data_list$outcome_data,
                      regimen_data=simulated_data_list$regimen_data,
                      baseline_data=simulated_data_list$sim_baseline_covariates,
                      timevar_data=simulated_data_list$sim_time_covariates,
                      det.Q.function=NULL,# now build-in
                      SL.library="glmnet",
                      SL.cvControl=list(selector="undersmooth",alpha=1),
                      verbose=TRUE)
    }),
    tar_target(table_elastic_estimate,{
        summary(elastic_estimate[[1]][[1]]$Ltmle_fit)
    })
)


######################################################################
### realistic_targets.R ends here
