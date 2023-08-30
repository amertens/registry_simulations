### get_simulated_data_list.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 17 2023 (14:14) 
## Version: 
## Last-Updated: Jul 17 2023 (15:23) 
##           By: Andrew Mertens
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log: Add cleaning of data so dementia occurs before death in 
### the same node, and if dementia/death occurs before
#----------------------------------------------------------------------
## 
### Code:
# add pnr number and then split wide format into a list of sources:
# treatment, outcome, covariates
get_simulated_data_list <- function(simulated_data, time_horizon=10){
    simulated_data[,pnr:=1:.N]
    names_time_covariates = c("heart.failure","renal.disease","chronic.pulmonary.disease","any.malignancy","ischemic.heart.disease","myocardial.infarction","hypertension","stroke","bb","ccb","rasi","thiazid","loop","mra","copd_med")
    sim_time_covariates <- simulated_data[,c("pnr",c(sapply(names_time_covariates,paste0,"_",0:(time_horizon-1)))),with = FALSE]
    sim_outcome<- simulated_data[,grep("pnr|dementia_|Censored|Dead", names(simulated_data)), with = FALSE]
    ## fix outcome to value 1 after first occurrence
    for (i in 1:time_horizon){
        for (j in ((i+1):time_horizon)){
            set(sim_outcome,i=which(sim_outcome[[paste0("dementia_",i)]]==1),j=paste0("dementia_",j),value=1)
        }
    }
    sim_baseline_covariates=simulated_data[,c("pnr","sex","agegroups","education","income","diabetes_duration"),with=FALSE]
    list(regimen_data = list("GLP1RA" = simulated_data[,grep("pnr|GLP1RA", names(simulated_data)), with = FALSE]),
         outcome_data = list(dementia=sim_outcome),
         sim_baseline_covariates = sim_baseline_covariates,
         sim_time_covariates = sim_time_covariates)
}


######################################################################
### get_simulated_data_list.R ends here
