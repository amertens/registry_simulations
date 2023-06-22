

#Set dementia after death to NA
clean_sim_data <- function(d, N_time=10){
  
  d<- data.table(d)
  
  for(i in 1:(N_time+1)){
    j=i+1
    d[is.na(get(paste0("event_dementia_",i))), (paste0("event_dementia_",j)):=NA]
    d[get(paste0("event_dementia_",i))==1, (paste0("event_dementia_",j)):=1]
    d[get(paste0("event_death_",i))==1, (paste0("event_death_",j)):=1]
    d[get(paste0("event_death_",i))==1, (paste0("event_dementia_",j)):=NA]
  }
  return(d)
}