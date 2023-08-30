



run_ltmle_sim <- function(sim_d_list, 
                          Ncores=4,
                          Niter=200,
                          time_horizon=3,
                          name_outcome = "event_dementia",
                          name_regimen = "glp1",
                          name_censoring = "censor",
                          censored_label = "censored",
                          name_comp.event = "event_death",
                          baseline_variables=baseline_vars,
                          longituninal_covariates=long_covariates,
                          treatment_vars="glp1",
                          outcome_vars=c("event_dementia","censor","event_death"),
                          SL.library = "glm",
                          SL.cvControl=NULL,
                          deterministic.Q.function = det.Q.function,
                          #Qint=F, # need to implement
                          #gcomp, # need to implement
                          test = FALSE,
                          #gbound = c(0.01, 1), #Need to implement
                          #varmethod = "ic", #Need to implement
                          concurrentY=TRUE,
                          Markov_vars=NULL){
  library(parallel)
  library(doParallel)
  registerDoParallel(cores=Ncores)
  
  set.seed(12345)
  resdf = NULL
  resdf <- foreach(i = 1:Niter, .combine = 'bind_rows', .errorhandling = 'remove') %dopar% {
    res <- fit <- NULL
    try(fit <- run_Ltmle(d=sim_d_list[[i]],
                         time_horizon=time_horizon,
                         name_outcome = name_outcome,
                         name_regimen = name_regimen,
                         name_censoring = name_censoring,
                         censored_label = censored_label,
                         name_comp.event = name_comp.event,
                         baseline_variables=baseline_vars,
                         longituninal_covariates=long_covariates,
                         treatment_vars=treatment_vars,
                         outcome_vars=outcome_vars,
                         SL.library = SL.library,
                         SL.cvControl = SL.cvControl,
                         deterministic.Q.function = deterministic.Q.function,
                         #Qint=F, # need to implement
                         #gcomp, # need to implement
                         test = test,
                         #gbound = c(0.01, 1), #Need to implement
                         #varmethod = "ic", #Need to implement
                         concurrentY=concurrentY,
                         Markov_vars=Markov_variables
                         ))
    try(res <- clean_ltmle_res(fit=fit, analysis_name="", iteration=i))
    return(res)
  }
  resdf
  return(resdf)
}

