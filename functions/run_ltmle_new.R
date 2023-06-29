run_ltmle_new <- function(d,
                          time_horizon=3,
                          baseline_vars=baseline_vars,
                          long_covariates=long_covariates,
                          treatment_vars="glp1",
                          outcome_vars=c("event_dementia","censor","event_death"),
                          N_time = 11, #number of time points you want to look at
                          SL.library = c("SL.glmnet"),
                          resdf=NULL,
                          Qint=F,
                          gcomp=F,
                          det.Q=T,
                          gbound = c(0.01, 1),
                          override_function=SuperLearner_override,
                          varmethod = "tmle", #variance method
                          label="",
                          glm=FALSE,
                          id=NULL){
  
  d_wide_list <- tar_read(sim_data)
  
  
    d <- data.table(d_wide_list[[1]])
    outcome_data <- d[,c("pnr",grep(paste0(outcome_vars, collapse ="|"),names(d), value = TRUE)), with = FALSE]
    treatment_data <- d[,c("pnr",grep(treatment_vars,names(d), value = TRUE)), with = FALSE]
    baseline_data <- d[,.(pnr,baseline_vars)]
    timevar_data <- d[,c("pnr",grep(paste0(long_covariates, collapse ="|"),names(d), value = TRUE)), with = FALSE] 
  
  
    pl <- prepare_Ltmle(outcome_data = list("event_dementia"=outcome_data),
                        regimen_data = treatment_data,
                        baseline_data = baseline_data,
                        timevar_data = timevar_data,
                        time_horizon = time_horizon,
                        deterministic.Q.function = det.Q.function,
                        name_outcome = "event_dementia",
                        name_regimen = "glp1",
                        name_censoring = "censor",
                        censored_label = "censored",
                        name_comp.event = "event_death",
                        Markov = NULL, #set to true?
                        subset_id = NULL,
                        SL.library="glm",
                        test = FALSE,
                        abar = list(rep(1,time_horizon),rep(0,time_horizon)))
  
  
    pl$verbose=1L 
    fit <- do.call(Ltmle, pl)
    ## fit <- do.call(ltmle::ltmle, pl)
    

  print(summary(fit))
  
  # [1] "IC"                "msm"               "beta"              "cum.g"             "cum.g.unbounded"   "fit"              
  # [7] "variance.estimate" "beta.iptw"         "IC.iptw"           "Qstar"             "cum.g.used"        "gcomp"            
  # [13] "formulas"          "binaryOutcome"     "transformOutcome"  "survivalOutcome"   "call"              "info" 
  
 
  return(fit)
}
