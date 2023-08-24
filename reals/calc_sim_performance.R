




calc_sim_performance <- function(res, truth, time=10, mean=TRUE){
  
  if(class(res)[1]=="list"){
    res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")
  }else{
    if(is.null(res$analysis)){
      res$analysis = ""
    }
  }
  
  if(mean){
    trueYa1 <- truth[time,2]
    trueRR <- truth[time,3]
    trueRD <- truth[time,4]
  }else{
    trueYa1 <- truth[time,5]
    trueRR <- truth[time,6]
    trueRD <- truth[time,7]   
  }

  

  res_Ya1 <- res %>% filter(Target_parameter=="Risk(A=1)") %>%
    group_by(Estimator , estimator) %>% 
    summarise(N_reps=n(),
              abs_bias_Ya1=mean(abs(estimate-trueYa1)),
              estimator_variance_Ya1=mean(((estimate)-mean((estimate)))^2),
              mean_variance_Ya1=mean((std.err)^2),
              bias_se_ratio_Ya1=abs_bias_Ya1/sqrt(mean_variance_Ya1),
              coverage_Ya1=mean(lower<=trueYa1 & trueYa1<=upper)*100,
              O_coverage_Ya1=mean(estimate-1.96*sd(estimate)< trueYa1 & trueYa1 < estimate+1.96*sd(estimate))*100
    )
    resRD <- res %>% filter(Target_parameter=="ATE") %>%
    group_by(Estimator, estimator) %>% 
    summarise(N_reps=n(),
              abs_bias_RD=mean(abs(estimate-trueRD)),
              estimator_variance_RD=mean(((estimate)-mean((estimate)))^2),
              mean_variance_RD=mean((std.err)^2),
              bias_se_ratio_RD=abs_bias_RD/sqrt(mean_variance_RD),
              coverage_RD=mean(lower<=trueRD & trueRD<=upper)*100,
              O_coverage_RD=mean(estimate-1.96*sd(estimate)< trueRD & trueRD < estimate+1.96*sd(estimate))*100
    )
  
  resRR <- res %>% filter(Target_parameter=="RelativeRisk") %>%
    group_by(Estimator, estimator) %>% 
    summarise(abs_log_bias_RR=mean(abs(log(estimate)-log(trueRR))),
              estimator_variance_RR=mean(((estimate)-mean((estimate)))^2),
              mean_variance_RR=mean((std.err)^2),
              bias_se_ratio_RR=abs_log_bias_RR/sqrt(mean_variance_RR),
              coverage_RR=mean(lower<=trueRR & trueRR<=upper)*100,
              O_coverage_RR=mean(log(estimate)-1.96*sd(log(estimate))< log(trueRR) & 
                                   log(trueRR) < log(estimate)+1.96*sd(log(estimate)))*100
    )
  
  tab = merge(res_Ya1, resRD, by=c("Estimator","estimator"))
  tab = merge(tab, resRR, by=c("Estimator","estimator"))
  
  return(tab)
}
