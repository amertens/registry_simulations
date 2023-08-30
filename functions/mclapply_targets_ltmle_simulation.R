

mclapply_targets_ltmle_simulation <- function(seeds, n_df=100000, n_cores=50,
                                              library="glm",
                                              SL.Control=NULL,
                                              time=10,
                                              gbounds=c(0.01,1),
                                              null_sim=FALSE,
                                              n_bootstrap_samples=0,
                                              Markov_variables=NULL,
                                              estimator="",
                                              tmle_var=FALSE){


  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("seeds","run_targets_ltmle_simulation",  "library",
                      "SL.Control","null_sim",
                      "n_bootstrap_samples",
                      "Markov_variables",
                      "gbounds",
                      "n_df",
                      "time","tmle_var"), envir=environment())
  
  clusterEvalQ(cl,lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
    do.call("require", list(X)) 
  }))


  res=parLapply(cl=cl,seeds, function(z) run_targets_ltmle_simulation(seed=z,
                                                            library=library,
                                                            null_sim=null_sim,
                                                            SL.Control=SL.Control,
                                                            n_bootstrap_samples=n_bootstrap_samples,
                                                            Markov_variables=Markov_variables,
                                                            gbounds=gbounds,
                                                            n=n_df,
                                                            time=time,
                                                            tmle_var=tmle_var))
  stopCluster(cl)

  drop <- rep(FALSE, length(res))
  for(i in 1:length(res)){
    if(class(res[[i]])[1]!="data.table"){
      drop[i] <- TRUE
    }
  }
  
  res[drop] <- NULL
  
  res=rbindlist(res, fill=TRUE)
  res$estimator=estimator
  return(res)
  print(Sys.time()-start) 
}
