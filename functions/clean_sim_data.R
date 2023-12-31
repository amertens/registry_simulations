

#Set dementia after death to NA
clean_sim_data <- function(d, N_time){
  
  d<- data.table(d)
  
  for(i in 1:(N_time+1)){
    j=i+1
    d[get(paste0("dementia_",i))==1, (paste0("dementia_",j)):=1]
    d[get(paste0("Dead_",i))==1, (paste0("Dead_",j)):=1]
  }
  
  dementia.nodes<- grep("dementia_",names(d))
  death.nodes<- grep("Dead_",names(d))
  d[, sum_death :=rowSums(.SD,na.rm=T), .SDcols = death.nodes]
  d[, sum_dementia :=rowSums(.SD,na.rm=T), .SDcols = dementia.nodes]
  d[sum_death > sum_dementia, (dementia.nodes) := replace(.SD, .SD == 1, 0), .SDcols = dementia.nodes]
  d[sum_death < sum_dementia, (death.nodes) := replace(.SD, .SD == 1, 0), .SDcols = death.nodes]
  d[sum_death== sum_dementia, (death.nodes) := replace(.SD, .SD == 1, 0), .SDcols = death.nodes]
  return(d)

}