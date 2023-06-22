

generate_data <- function(cc, seed=12345, reps=sim_reps, n=dataset_N, N_time=10, clean_death=FALSE){
  
  u <- synthesizeDD(cc)

  
  set.seed(seed)
  sim_list <- NULL

  for(i in 1:reps){
    cat("\ni: ",i,"\n")
    
    d <- sim(u,n)
    
    if(clean_death){
      #set dementia after death to NA
      d<- clean_sim_data(d, N_time = 10)      
    }
    
    #add identifiers
    d$pnr <- 1:nrow(d)
    
    sim_list[[i]] <- d
    gc()
  }
  
  return(sim_list)
}