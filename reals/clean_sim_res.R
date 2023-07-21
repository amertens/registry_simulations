

#Set dementia after death to NA
clean_sim_res <- function(res){
  res_tab=NULL
  for(i in 1:(length(res)/4)){
    temp = summary(res[[i + (3*(i-1))]]$GLP1RA$Ltmle_fit)
    res_tab <- bind_rows(res_tab, temp)
  }
  
  return(res_tab)
}