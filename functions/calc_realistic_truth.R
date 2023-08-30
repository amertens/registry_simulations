

# cc <- coefs <- tar_read(coefs)
# time_horizon=10
# nsamp=1000000


calc_realistic_truth <- function(A_name = "glp1", nsamp=100000, return_data=FALSE){
  
  source("data/coefs.txt")
  modelA1 = get_lava_model(time_horizon = 10, coefs = coefs, counterfactual_A=1)
  modelA0 = get_lava_model(time_horizon = 10, coefs = coefs, counterfactual_A=0)
  d.always = simulate_data(lava_model = modelA1, n = nsamp)
  d.never = simulate_data(lava_model = modelA0, n = nsamp)

  d.always <- clean_sim_data(d.always, N_time=10)
  d.never <- clean_sim_data(d.never, N_time=10)

  Ya1_t1 <- mean(d.always$dementia_1,na.rm=T)
  Ya1_t2 <- mean(d.always$dementia_2,na.rm=T)
  Ya1_t3 <- mean(d.always$dementia_3,na.rm=T)
  Ya1_t4 <- mean(d.always$dementia_4,na.rm=T)
  Ya1_t5 <- mean(d.always$dementia_5,na.rm=T)
  Ya1_t6 <- mean(d.always$dementia_6,na.rm=T)
  Ya1_t7 <- mean(d.always$dementia_7,na.rm=T)
  Ya1_t8 <- mean(d.always$dementia_8,na.rm=T)
  Ya1_t9 <- mean(d.always$dementia_9,na.rm=T)
  Ya1_t10 <- mean(d.always$dementia_10,na.rm=T)
  
    
  tRR1 <- mean(d.always$dementia_1,na.rm=T)/mean(d.never$dementia_1,na.rm=T)
  tRR2 <- mean(d.always$dementia_2,na.rm=T)/mean(d.never$dementia_2,na.rm=T)
  tRR3 <- mean(d.always$dementia_3,na.rm=T)/mean(d.never$dementia_3,na.rm=T)
  tRR4 <- mean(d.always$dementia_4,na.rm=T)/mean(d.never$dementia_4,na.rm=T)
  tRR5 <- mean(d.always$dementia_5,na.rm=T)/mean(d.never$dementia_5,na.rm=T)
  tRR6 <- mean(d.always$dementia_6,na.rm=T)/mean(d.never$dementia_6,na.rm=T)
  tRR7 <- mean(d.always$dementia_7,na.rm=T)/mean(d.never$dementia_7,na.rm=T)
  tRR8 <- mean(d.always$dementia_8,na.rm=T)/mean(d.never$dementia_8,na.rm=T)
  tRR9 <- mean(d.always$dementia_9,na.rm=T)/mean(d.never$dementia_9,na.rm=T)
  tRR10 <- mean(d.always$dementia_10,na.rm=T)/mean(d.never$dementia_10,na.rm=T)
  
  tRD1 <- mean(d.always$dementia_1,na.rm=T) - mean(d.never$dementia_1,na.rm=T)
  tRD2 <- mean(d.always$dementia_2,na.rm=T) - mean(d.never$dementia_2,na.rm=T)
  tRD3 <- mean(d.always$dementia_3,na.rm=T) - mean(d.never$dementia_3,na.rm=T)
  tRD4 <- mean(d.always$dementia_4,na.rm=T) - mean(d.never$dementia_4,na.rm=T)
  tRD5 <- mean(d.always$dementia_5,na.rm=T) - mean(d.never$dementia_5,na.rm=T)
  tRD6 <- mean(d.always$dementia_6,na.rm=T) - mean(d.never$dementia_6,na.rm=T)
  tRD7 <- mean(d.always$dementia_7,na.rm=T) - mean(d.never$dementia_7,na.rm=T)
  tRD8 <- mean(d.always$dementia_8,na.rm=T) - mean(d.never$dementia_8,na.rm=T)
  tRD9 <- mean(d.always$dementia_9,na.rm=T) - mean(d.never$dementia_9,na.rm=T)
  tRD10 <- mean(d.always$dementia_10,na.rm=T) - mean(d.never$dementia_10,na.rm=T)
  
  truth_df <- data.frame(time=1:10, 
                         Ya1=c(  Ya1_t1,
                                 Ya1_t2,
                                 Ya1_t3,
                                 Ya1_t4,
                                 Ya1_t5,
                                 Ya1_t6,
                                 Ya1_t7,
                                 Ya1_t8,
                                 Ya1_t9,
                                 Ya1_t10),
                         RR=c(tRR1,tRR2,tRR3,tRR4,tRR5,tRR6,tRR7,tRR8,tRR9,tRR10), 
                         RD=c(tRD1,tRD2,tRD3,tRD4,tRD5,tRD6,tRD7,tRD8,tRD9,tRD10))

  if(return_data){
    return(list(truth_df=truth_df, dA1=d.always, dA0=d.never))
  }else{
    return(truth_df)
  }
}


average_truth <- function(truth){
  truth <- truth %>% group_by(time) %>%
    summarise(meanYa1=mean(Ya1), meanRR=mean(RR), meanRD=mean(RD),
              medianYa1=median(Ya1), medianRR=median(RR), medianRD=median(RD)) %>%
    as.data.frame()
  return(truth)
}
