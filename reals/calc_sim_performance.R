

calc_sim_performance <- function(res, truth, time=10){
  
  trueRR <- truth[time,2]
  trueRD <- truth[time,3]
  
  resRD <- res$estimate[res$Target_parameter=="ATE"]
  resRR <- res$estimate[res$Target_parameter=="RelativeRisk"]
  
  
  O_coverage_RD = mean(resRD-1.96*sd(resRD)< trueRD & trueRD < resRD+1.96*sd(resRD))  
  O_coverage_RR = mean(resRR-1.96*sd(resRR)< (trueRR) & trueRR < resRR+1.96*sd(resRR))  
  
  tab = data.frame(O_coverage_RD=O_coverage_RD,  O_coverage_RR=O_coverage_RR)
  
  return(tab)
}


