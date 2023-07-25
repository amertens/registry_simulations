




calc_sim_performance <- function(res, truth, time=10){
  
  if(class(res)=="list"){
    res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")
  }else{
    res$analysis = ""
  }
  
  trueRR <- truth[time,2]
  trueRD <- truth[time,3]
  
  resRD <- res %>% filter(Target_parameter=="ATE") %>%
    group_by(analysis) %>% 
    summarise(abs_bias_RD=mean(abs(estimate-trueRD)),
              variance_RD=mean(((estimate)-mean((estimate)))^2),
              coverage_RD=mean(lower<=trueRD & trueRD<=upper)*100,
              O_coverage_RD=mean(estimate-1.96*sd(estimate)< trueRD & trueRD < estimate+1.96*sd(estimate))*100
    )
  
  resRR <- res %>% filter(Target_parameter=="RelativeRisk") %>%
    group_by(analysis) %>% 
    summarise(abs_bias_RR=mean(abs(log(estimate)-log(trueRR))),
              variance_RR=mean(((estimate)-mean((estimate)))^2),
              coverage_RR=mean(lower<=trueRR & trueRR<=upper)*100,
              O_coverage_RR=mean(log(estimate)-1.96*sd(log(estimate))< log(trueRR) & 
                                   log(trueRR) < log(estimate)+1.96*sd(log(estimate)))*100
    )
  
  tab = merge(resRD, resRR, by=c("analysis"))
  
  return(tab)
}
