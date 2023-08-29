

cc <- coefs <- tar_read(coefs)
time_horizon=10
nsamp=1000000


calc_truth <- function(cc, A_name = "glp1", seed=12345,  nsamp=100000){
  
 #To do: check if seed is a vector and then duplicate and average the calculated results
  #if(lenght(seed)>1){
  
  #}else{
  
  set.seed(seed)
  u.always <- synthesizeDD(cc, A=1)
  d.always  <- sim(u.always, nsamp)
  
  set.seed(seed)
  u.never <- synthesizeDD(cc, A=0)
  d.never <- sim(u.never, nsamp)
  
  d.always<- clean_sim_data(d.always, N_time=10)
  d.never<- clean_sim_data(d.never, N_time=10)
  
  #to do: update so a flexible number of time point truths can be calculated
  
  tRR1 <- mean(d.always$event_dementia_1,na.rm=T)/mean(d.never$event_dementia_1,na.rm=T)
  tRR2 <- mean(d.always$event_dementia_2,na.rm=T)/mean(d.never$event_dementia_2,na.rm=T)
  tRR3 <- mean(d.always$event_dementia_3,na.rm=T)/mean(d.never$event_dementia_3,na.rm=T)
  tRR4 <- mean(d.always$event_dementia_4,na.rm=T)/mean(d.never$event_dementia_4,na.rm=T)
  tRR5 <- mean(d.always$event_dementia_5,na.rm=T)/mean(d.never$event_dementia_5,na.rm=T)
  tRR6 <- mean(d.always$event_dementia_6,na.rm=T)/mean(d.never$event_dementia_6,na.rm=T)
  tRR7 <- mean(d.always$event_dementia_7,na.rm=T)/mean(d.never$event_dementia_7,na.rm=T)
  tRR8 <- mean(d.always$event_dementia_8,na.rm=T)/mean(d.never$event_dementia_8,na.rm=T)
  tRR9 <- mean(d.always$event_dementia_9,na.rm=T)/mean(d.never$event_dementia_9,na.rm=T)
  tRR10 <- mean(d.always$event_dementia_10,na.rm=T)/mean(d.never$event_dementia_10,na.rm=T)
  
  tRD1 <- mean(d.always$event_dementia_1,na.rm=T) - mean(d.never$event_dementia_1,na.rm=T)
  tRD2 <- mean(d.always$event_dementia_2,na.rm=T) - mean(d.never$event_dementia_2,na.rm=T)
  tRD3 <- mean(d.always$event_dementia_3,na.rm=T) - mean(d.never$event_dementia_3,na.rm=T)
  tRD4 <- mean(d.always$event_dementia_4,na.rm=T) - mean(d.never$event_dementia_4,na.rm=T)
  tRD5 <- mean(d.always$event_dementia_5,na.rm=T) - mean(d.never$event_dementia_5,na.rm=T)
  tRD6 <- mean(d.always$event_dementia_6,na.rm=T) - mean(d.never$event_dementia_6,na.rm=T)
  tRD7 <- mean(d.always$event_dementia_7,na.rm=T) - mean(d.never$event_dementia_7,na.rm=T)
  tRD8 <- mean(d.always$event_dementia_8,na.rm=T) - mean(d.never$event_dementia_8,na.rm=T)
  tRD9 <- mean(d.always$event_dementia_9,na.rm=T) - mean(d.never$event_dementia_9,na.rm=T)
  tRD10 <- mean(d.always$event_dementia_10,na.rm=T) - mean(d.never$event_dementia_10,na.rm=T)
  
  truth_df <- data.frame(time=1:10, RR=c(tRR1,tRR2,tRR3,tRR4,tRR5,tRR6,tRR7,tRR8,tRR9,tRR10), RD=c(tRD1,tRD2,tRD3,tRD4,tRD5,tRD6,tRD7,tRD8,tRD9,tRD10))
 #}
  return(list(truth_df=truth_df, dA1=d.always, dA0=d.never))
}