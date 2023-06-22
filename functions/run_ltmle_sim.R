

run_ltmle_sim <- function(){
  
  d_wide_list <- tar_read(sim_data)
  
  library(parallel)
  library(doParallel)
  registerDoParallel(cores=4)
  
  
  resdf_DetQ_ic <- foreach(i = 1:length(d_wide_list), .combine = 'bind_rows', .errorhandling = 'remove') %dopar% {
    res <- NULL
    try(res <- run_ltmle(d_wide_list[[i]], resdf=NULL, Qint=FALSE, det.Q =FALSE, varmethod = "ic"))
    return(res)
  }
  
  
}

