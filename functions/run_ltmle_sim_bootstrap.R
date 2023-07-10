

# {
# sim_d_list <- tar_read(sim_data)
# Ncores=4
# Niter=2
# Nbootstrap=2
# time_horizon=3
# name_outcome = "event_dementia"
# name_regimen = "glp1"
# name_censoring = "censor"
# censored_label = "censored"
# name_comp.event = "event_death"
# baseline_variables=baseline_vars
# longituninal_covariates=long_covariates
# treatment_vars="glp1"
# outcome_vars=c("event_dementia","censor","event_death")
# SL.library = "glm"
# SL.cvControl=NULL
# deterministic.Q.function = det.Q.function
# test = FALSE
# Markov_vars=Markov_variables
# }

run_ltmle_sim_bootstrap <- function(sim_d_list, 
                          Ncores=4,
                          Niter=200,
                          Nbootstrap=200,
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
  
  i<-j<-1
  resdf_boot = NULL
  
  for(i in 1:Nbootstrap){
    
    cat(i,"\n")
    res_df <- NULL
    res_df <- foreach(j = 1:Niter, .combine = 'bind_rows', .errorhandling = 'remove') %dopar% {
      
      set.seed(j)
      d <- sim_d_list[[j]]
      d$id <- 1:nrow(d)
      dboot <- d[sample(.N, nrow(d),replace=TRUE)]
      #hack: remake unique PNR's for merge
      dboot$pnr <- 1:nrow(dboot)
      res <- fit <- NULL
      
      #TO DO
      #check passing of ID arg
      try(fit <- run_Ltmle(d= dboot,
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
                           id=dboot$id,
                           concurrentY=concurrentY,
                           Markov_vars=NULL
      ))
      try(res <- clean_ltmle_res(fit=fit, analysis_name="", iteration=i))
      return(res)
      
    }
    res_df
    
    gc()
    res_df$iteration <- i
    resdf_boot <- bind_rows(resdf_boot, res_df)
  }
  
  return(resdf_boot)
}

