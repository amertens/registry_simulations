
rm(list=ls())
# library(tidyverse)
 # simulated_data = readRDS(paste0(here::here(),"/data/sim_data/alt/simulated_data_",1,".rds"))
 # simulated_data_list=get_simulated_data_list(simulated_data)
 # ntime=10
# colnames(simulated_data)
# simulated_data <- simulated_data[1:1000,]
# df <- simulated_data %>% dplyr::select(sex, hypertension_0, GLP1RA_0, dementia_1, hypertension_1,  GLP1RA_1, dementia_2, Censored_1, Censored_2, Dead_1, Dead_2)
# dput(df)

source(here::here("reproducible_example_functions.R"))
library(data.table)

#------------------------------------------------------------------------------
# prep data
#------------------------------------------------------------------------------

head(d)
simulated_data=data.table(d)
simulated_data[,pnr:=1:.N]

names_time_covariates = c("hypertension")
ntime=2

sim_time_covariates <- simulated_data[,c("pnr",c(sapply(names_time_covariates,paste0,"_",0:(ntime-1)))),with = FALSE]
sim_outcome<- simulated_data[,grep("pnr|dementia_|Censored|Dead", names(simulated_data)), with = FALSE]
## fix outcome to value 1 after first occurrence
for (i in 1:ntime){
  for (j in ((i+1):ntime)){
    set(sim_outcome,i=which(sim_outcome[[paste0("dementia_",i)]]==1),j=paste0("dementia_",j),value=1)
  }
}
sim_baseline_covariates=simulated_data[,c("pnr","sex"),with=FALSE]
simulated_data_list=list(regimen_data = list("GLP1RA" = simulated_data[,grep("pnr|GLP1RA", names(simulated_data)), with = FALSE]),
     outcome_data = list(dementia=sim_outcome),
     sim_baseline_covariates = sim_baseline_covariates,
     sim_time_covariates = sim_time_covariates)





#------------------------------------------------------------------------------
# Function with bug called from CalcIPTW, FitPooledMSM, UpdateQ (via FixedTimeTMLE) with Q* as outcome
#------------------------------------------------------------------------------

ltmle.glm<- function (formula, family, data, weights){
  
  #Old code: length(weights)==NROW(data) means speedglm without weights is used for all estimates with this function
  #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  if (length(weights)==NROW(data)||is.null(weights)) {
  #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    try.result <- try(m <- speedglm::speedglm(formula = formula,
                                              family = family,
                                              data = data,
                                              maxit = 100),
                      silent = TRUE)
    if (inherits(try.result, "try-error")) {
      ShowGlmMessage()
      m <- glm(formula = formula,
               family = family,
               data = data,
               control = glm.control(maxit = 100))
    }
  }
  else {
    m <- glm(formula = formula,
             family = family,
             data = data.frame(data,weights),
             weights = weights,
             control = glm.control(maxit = 100))
  }
  return(m)
}

#------------------------------------------------------------------------------
# Run comparison of glm and glmnet with bug
#------------------------------------------------------------------------------


set.seed(1234)
res_glm_old = run_ltmle(name_outcome="dementia",
                time_horizon=ntime,
                test=FALSE,
                outcome_data=simulated_data_list$outcome_data,
                regimen_data=simulated_data_list$regimen_data,
                baseline_data=simulated_data_list$sim_baseline_covariates,
                timevar_data=simulated_data_list$sim_time_covariates,
                det.Q.function=NULL,# now build-in
                gbounds=c(0.1,1),
                B_bootstrap_samples=0,
                SL.library="glm",
                Markov="hypertension",
                SL.cvControl=NULL,
                tmle_var=FALSE,
                verbose=TRUE)

set.seed(1234)
res_glmnet_old = run_ltmle(name_outcome="dementia",
                        time_horizon=ntime,
                        test=FALSE,
                        outcome_data=simulated_data_list$outcome_data,
                        regimen_data=simulated_data_list$regimen_data,
                        baseline_data=simulated_data_list$sim_baseline_covariates,
                        timevar_data=simulated_data_list$sim_time_covariates,
                        det.Q.function=NULL,# now build-in
                        gbounds=c(0.1,1),
                        B_bootstrap_samples=0,
                        SL.library="glmnet",
                        Markov="hypertension",
                        SL.cvControl=list(selector="undersmooth",alpha=1),
                        tmle_var=FALSE,
                        verbose=TRUE)

# #tmle results vary between glmnet and glm-fit ltmle
# summary(res_glm_old$time_horizon_2$GLP1RA$Ltmle_fit)
# summary(res_glmnet_old$time_horizon_2$GLP1RA$Ltmle_fit)
# 
# 
# #tmle results vary between glmnet and glm-fit ltmle
# summary(res_glm_old$time_horizon_2$GLP1RA$Ltmle_fit, estimator = "iptw")
# summary(res_glmnet_old$time_horizon_2$GLP1RA$Ltmle_fit, estimator = "iptw")

#tmle results vary between glmnet and glm-fit ltmle
summary(res_glm_old$time_horizon_10$GLP1RA$Ltmle_fit)
summary(res_glmnet_old$time_horizon_10$GLP1RA$Ltmle_fit)


#tmle results vary between glmnet and glm-fit ltmle
summary(res_glm_old$time_horizon_10$GLP1RA$Ltmle_fit, estimator = "iptw")
summary(res_glmnet_old$time_horizon_10$GLP1RA$Ltmle_fit, estimator = "iptw")

#------------------------------------------------------------------------------
# FIXED Function called from CalcIPTW, FitPooledMSM, UpdateQ (via FixedTimeTMLE) with Q* as outcome
#------------------------------------------------------------------------------

ltmle.glm<- function (formula, family, data, weights){
  
  #new code (also fixed in LTMLE for breakfast)
  #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  if (length(weights)!=NROW(data)||is.null(weights)) {
    #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    try.result <- try(m <- speedglm::speedglm(formula = formula,
                                              family = family,
                                              data = data,
                                              maxit = 100),
                      silent = TRUE)
    if (inherits(try.result, "try-error")) {
      ShowGlmMessage()
      m <- glm(formula = formula,
               family = family,
               data = data,
               control = glm.control(maxit = 100))
    }
  }
  else {
    m <- glm(formula = formula,
             family = family,
             data = data.frame(data,weights),
             weights = weights,
             control = glm.control(maxit = 100))
  }
  return(m)
}


set.seed(1234)
res_glm_new = run_ltmle(name_outcome="dementia",
                        time_horizon=ntime,
                        test=FALSE,
                        outcome_data=simulated_data_list$outcome_data,
                        regimen_data=simulated_data_list$regimen_data,
                        baseline_data=simulated_data_list$sim_baseline_covariates,
                        timevar_data=simulated_data_list$sim_time_covariates,
                        det.Q.function=NULL,# now build-in
                        gbounds=c(0.1,1),
                        B_bootstrap_samples=0,
                        SL.library="glm",
                        Markov="hypertension",
                        SL.cvControl=NULL,
                        tmle_var=FALSE,
                        verbose=TRUE)

set.seed(1234)
res_glmnet_new = run_ltmle(name_outcome="dementia",
                           time_horizon=ntime,
                           test=FALSE,
                           outcome_data=simulated_data_list$outcome_data,
                           regimen_data=simulated_data_list$regimen_data,
                           baseline_data=simulated_data_list$sim_baseline_covariates,
                           timevar_data=simulated_data_list$sim_time_covariates,
                           det.Q.function=NULL,# now build-in
                           gbounds=c(0.1,1),
                           B_bootstrap_samples=0,
                           SL.library="glmnet",
                           Markov="hypertension",
                           SL.cvControl=list(selector="undersmooth",alpha=1),
                           tmle_var=FALSE,
                           verbose=TRUE)

#tmle results vary between glmnet and glm-fit ltmle
summary(res_glm_new$time_horizon_2$GLP1RA$Ltmle_fit)
summary(res_glmnet_new$time_horizon_2$GLP1RA$Ltmle_fit)


#tmle results vary between glmnet and glm-fit ltmle
summary(res_glm_new$time_horizon_2$GLP1RA$Ltmle_fit, estimator = "iptw")
summary(res_glmnet_new$time_horizon_2$GLP1RA$Ltmle_fit, estimator = "iptw")


