
# seeds=c(1,2)
# n_df=10000
# n_cores=5
# library="glm"
# SL.Control=NULL
# time=10
# gbounds=c(0.01,1)
# n_bootstrap_samples=0
# Markov_variables=NULL
# 
# estimator="undersmoothed ridge"
# seeds=seeds1
# gbounds=c(0,1)
# library="glmnet"
# SL.Control=list(selector="undersmooth",alpha=0)



mclapply_targets_ltmle_simulation <- function(seeds, n_df=100000, n_cores=50,
                                              library="glm",
                                              SL.Control=NULL,
                                              time=10,
                                              gbounds=c(0.01,1),
                                              n_bootstrap_samples=0,
                                              Markov_variables=NULL,
                                              estimator=""){
  
  start=Sys.time()
  # res= mclapply(seeds, function(z) run_targets_ltmle_simulation(seed=z,
  #                                                               library=library,
  #                                                               SL.Control=SL.Control,
  #                                                               n_bootstrap_samples=n_bootstrap_samples,
  #                                                               Markov_variables=Markov_variables,
  #                                                               gbounds=gbounds,
  #                                                               n=n_df,
  #                                                               time=time),
  #                                   mc.set.seed =FALSE,
  #                                   mc.cores=n_cores)

  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("seeds","run_targets_ltmle_simulation",  "library",
                      "SL.Control",
                      "n_bootstrap_samples",
                      "Markov_variables",
                      "gbounds",
                      "n_df",
                      "time"))
  clusterEvalQ(cl,lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
    do.call("require", list(X)) 
  }))

    #library(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","Publish","matrixStats","speedglm","parallel","caret","foreach","clustermq")))
  

  # test=run_targets_ltmle_simulation(seed=seeds[3],
  #                                   library=library,
  #                                   SL.Control=SL.Control,
  #                                   n_bootstrap_samples=n_bootstrap_samples,
  #                                   Markov_variables=Markov_variables,
  #                                   gbounds=gbounds,
  #                                   n=n_df,
  #                                   time=time)
  
  res=parLapply(cl=cl,seeds, function(z) run_targets_ltmle_simulation(seed=z,
                                                            library=library,
                                                            SL.Control=SL.Control,
                                                            n_bootstrap_samples=n_bootstrap_samples,
                                                            Markov_variables=Markov_variables,
                                                            gbounds=gbounds,
                                                            n=n_df,
                                                            time=time))
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
