




calc_sim_performance <- function(res, truth, time=10){
  
  if(class(res)[1]=="list"){
    res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")
  }else{
    if(is.null(res$analysis)){
      res$analysis = ""
    }
  }
  
  trueRR <- truth[time,2]
  trueRD <- truth[time,3]
  
  resRD <- res %>% filter(Target_parameter=="ATE") %>%
    group_by(estimator) %>% 
    summarise(N_reps=n(),
              abs_bias_RD=mean(abs(estimate-trueRD)),
              estimator_variance_RD=mean(((estimate)-mean((estimate)))^2),
              mean_variance_RD=mean((std.err)^2),
              bias_se_ratio=abs_bias_RD/sqrt(mean_variance_RD),
              coverage_RD=mean(lower<=trueRD & trueRD<=upper)*100,
              O_coverage_RD=mean(estimate-1.96*sd(estimate)< trueRD & trueRD < estimate+1.96*sd(estimate))*100
    )
  
  resRR <- res %>% filter(Target_parameter=="RelativeRisk") %>%
    group_by(estimator) %>% 
    summarise(abs_log_bias_RR=mean(abs(log(estimate)-log(trueRR))),
              estimator_variance_RR=mean(((estimate)-mean((estimate)))^2),
              mean_variance_RR=mean((std.err)^2),
              bias_se_ratio=abs_log_bias_RR/sqrt(mean_variance_RR),
              coverage_RR=mean(lower<=trueRR & trueRR<=upper)*100,
              O_coverage_RR=mean(log(estimate)-1.96*sd(log(estimate))< log(trueRR) & 
                                   log(trueRR) < log(estimate)+1.96*sd(log(estimate)))*100
    )
  
  tab = merge(resRD, resRR, by=c("estimator"))
  
  return(tab)
}
