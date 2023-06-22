
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
  print(fit)
  names(fit)
  summary(fit) 
  
  warn = getOption("warn")
  options(warn=-1)
  
  #clean competing events
  d <-clean_sim_data(d, N_time=N_time)
  
  if(!is.null(id)){
    baseline_vars <- c(baseline_vars,"id")
  }
  
  #Use only first N time points
  d <- d %>%
    dplyr::select(!!(baseline_vars),matches(paste0("_(",paste0(0:(N_time-1),collapse="|"),")$")))
  
  
  spec_ltmle <- spec_analysis_sim(data=d, c(long_covariates,"event_death_"),
                                  baseline_vars, N_time,
                                  Avars=c("glp1_"),
                                  Yvars=c("event_dementia_"),
                                  Cvars=c("censor_"))
  abar_spec = list(rep(1,N_time-1),rep(0,N_time-1))
  spec_ltmle$data <- spec_ltmle$data %>% select(ie_type,age_base,sex,code5txt,quartile_income,id, everything() )
  
  
  set.seed(12345)
  fit = NULL
  
  
  if(Qint){
    
    if(N_time==11){
      
      # qform = c(
      #   insulin_0="Q.kplus1 ~ 1",
      #   event_dementia_1="Q.kplus1 ~ 1",
      #   event_dementia_2="Q.kplus1 ~ 1",
      #   event_dementia_3="Q.kplus1 ~ 1",
      #   event_dementia_4="Q.kplus1 ~ 1",
      #   event_dementia_5="Q.kplus1 ~ 1",
      #   event_dementia_6="Q.kplus1 ~ 1",
      #   event_dementia_7="Q.kplus1 ~ 1",
      #   event_dementia_8="Q.kplus1 ~ 1",
      #   event_dementia_9="Q.kplus1 ~ 1",
      #   event_dementia_10="Q.kplus1 ~ 1"
      # )
      qform = c(
        insulin_1="Q.kplus1 ~ 1",
        insulin_2="Q.kplus1 ~ 1",
        insulin_3="Q.kplus1 ~ 1",
        insulin_4="Q.kplus1 ~ 1",
        insulin_5="Q.kplus1 ~ 1",
        insulin_6="Q.kplus1 ~ 1",
        insulin_7="Q.kplus1 ~ 1",
        insulin_8="Q.kplus1 ~ 1",
        insulin_9="Q.kplus1 ~ 1",
        insulin_10="Q.kplus1 ~ 1"
      )
    }
    if(N_time==4){
      qform = c(
        insulin_0="Q.kplus1 ~ 1",
        event_dementia_1="Q.kplus1 ~ 1",
        event_dementia_2="Q.kplus1 ~ 1",
        event_dementia_3="Q.kplus1 ~ 1")
    }
    if(N_time==2){
      qform = c(
        insulin_0="Q.kplus1 ~ 1",
        event_dementia_1="Q.kplus1 ~ 1")
    }
  }else{
    qform=NULL
  }
  
  
  if(det.Q){
    det.q.fun = det.Q.function
  }else{
    det.q.fun = NULL
  }
  
  if(!is.null(id)){
    id <- spec_ltmle$data[["id"]]
  }
  
  if(glm){
    
    try(fit <- ltmle(data=spec_ltmle$data,
                     Anodes = spec_ltmle$Anodes,
                     Cnodes = spec_ltmle$Cnodes,
                     Lnodes = spec_ltmle$Lnodes,
                     Ynodes = spec_ltmle$Ynodes,
                     gbound=gbound,
                     survivalOutcome = T,
                     abar = abar_spec,
                     gcomp=gcomp,
                     Qform = qform,
                     estimate.time=F,
                     deterministic.Q.function = det.q.fun,
                     SL.library = "glm",
                     variance.method = varmethod,
                     id=id
    ))
    
  }else{
    package_stub("SuperLearner", "SuperLearner", override_function, {
      testthatsomemore::package_stub("ltmle", "Estimate", Estimate_override, {
        try(fit <- ltmle(data=spec_ltmle$data,
                         Anodes = spec_ltmle$Anodes[spec_ltmle$Anodes!="glp1_0"],
                         Cnodes = spec_ltmle$Cnodes[spec_ltmle$Cnodes!="censor_0"],
                         Lnodes = spec_ltmle$Lnodes[spec_ltmle$Lnodes!="event_death_0"],
                         Ynodes = spec_ltmle$Ynodes[spec_ltmle$Ynodes!="event_dementia_0"],
                         gbound=gbound,
                         survivalOutcome = T,
                         abar = abar_spec,
                         gcomp=gcomp,
                         Qform = qform,
                         estimate.time=F,
                         deterministic.Q.function = det.q.fun,
                         SL.library = SL.library,
                         variance.method = varmethod,
                         id=id
        ))
      })})
    
  }
  
  
  if(!is.null(fit)){
    res <- summary(fit)
    res.iptw <- summary(fit, estimator="iptw")
    res.RR <- as.data.frame(res$effect.measures$RR)
    res.ate <- as.data.frame(res$effect.measures$ATE) %>% rename(ate.long.name=long.name,ate=estimate, ate.sd=std.dev , ate.pval=pvalue, ate.ci.lb=CI.2.5., ate.ci.ub=  CI.97.5., ate.log.std.err=log.std.err)
    
    res.RR.iptw <- as.data.frame(res.iptw$effect.measures$RR) %>% rename(iptw.long.name=long.name, iptw.estimate=estimate, iptw.sd=std.dev , iptw.pval=pvalue, iptw.ci.lb=CI.2.5., iptw.ci.ub=  CI.97.5., iptw.log.std.err=log.std.err)
    res.ate.iptw <- as.data.frame(res.iptw$effect.measures$ATE) %>% rename(iptw.ate.long.name=long.name, iptw.ate=estimate, iptw.ate.sd=std.dev , iptw.ate.pval=pvalue, iptw.ate.ci.lb=CI.2.5., iptw.ate.ci.ub=  CI.97.5., iptw.ate.log.std.err=log.std.err)
    
    res <- cbind(res.RR, res.ate, res.RR.iptw, res.ate.iptw)
    res$label <- label
  }
  if(!is.null(resdf)){
    res <- bind_rows(resdf, res)
  }
  
  options(warn=warn)
  return(res)
}