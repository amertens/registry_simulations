
run_Ltmle <- function(d,
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
                          concurrentY=FALSE, #Whether Yt is predicted from Lt (TRUE, used for sim data from coefficients) or Lt-1 (FALSE)
                          #Qint=F, # need to implement
                          #gcomp, # need to implement
                          test = FALSE,
                          #gbound = c(0.01, 1), #Need to implement
                          #varmethod = "ic", #Need to implement
                          id=NULL,
                          Markov_vars=Markov_variables){
  

    outcome_data <- d[,c("pnr",grep(paste0(outcome_vars, collapse ="|"),names(d), value = TRUE)), with = FALSE]
    treatment_data <- d[,c("pnr",grep(treatment_vars,names(d), value = TRUE)), with = FALSE]
    baseline_data <- d[,c("pnr",grep(paste0(baseline_variables, collapse ="|"),names(d), value = TRUE)), with = FALSE] 
    timevar_data <- d[,c("pnr",grep(paste0(longituninal_covariates, collapse ="|"),names(d), value = TRUE)), with = FALSE] 
  
    #browser()
    pl <- prepare_Ltmle_sim(outcome_data = list("event_dementia"=outcome_data),
                        regimen_data = treatment_data,
                        baseline_data = baseline_data,
                        timevar_data = timevar_data,
                        time_horizon = time_horizon,
                        deterministic.Q.function = det.Q.function,
                        name_outcome = name_outcome,
                        name_regimen = name_regimen,
                        name_censoring = name_censoring,
                        censored_label = censored_label,
                        name_comp.event = name_comp.event,
                        Markov = Markov_vars, #names of time-varying variables assumed to be following the markov process
                        concurrentY=concurrentY,
                        subset_id = NULL,
                        SL.library=SL.library, 
                        test = test,
                        abar = list(rep(1,time_horizon+ (1 * concurrentY)),rep(0,time_horizon+ (1 * concurrentY))))
  
    
    pl$verbose=1L 
    pl$id = NULL
    if (length(SL.cvControl)>0){pl$SL.cvControl <- SL.cvControl}
    # browser()
    # head(pl$data)
    
    fit <- do.call(estimate_Ltmle, pl)
    
    #NOTE:
    # use code from the run_ltmle function and the "keep" arguement to save regression coef and g weights

  return(fit)
}
