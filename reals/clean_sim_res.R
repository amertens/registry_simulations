

#Set dementia after death to NA
clean_sim_res <- function(res){
  res_tab=NULL
  res=res[grepl("Ltmle_fit",names(res))]
  for(i in 1:length(res)){
    temp = summary(res[[i]])
    temp$iteration = i
    res_tab <- bind_rows(res_tab, temp)
  }
  
  return(res_tab)
}