
rm(list=ls())
load(paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))

gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)

# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
#set
#Longitudinal variables possibly following the markov process
Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                   "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                   "bb","ccb","rasi","thiazid",
                   "loop","mra","copd_med"  )

 # n_cores=50
 # estimator="glm"
 # dataset_nums=1
 # library="glm"
 # n_df=100000
 # SL.Control=NULL
 # time=10
 # gbounds=c(0.01,1)
 # null_sim=FALSE
 # n_bootstrap_samples=0
 # Markov_variables=NULL
 # tmle_var=FALSE
 # library="glm"

mclapply_targets_ltmle_simulation2 <- function(dataset_nums, n_df=100000, n_cores=50,
         library="glm",
         SL.Control=NULL,
         time=10,
         gbounds=c(0.01,1),
         null_sim=FALSE,
         n_bootstrap_samples=0,
         Markov_variables=NULL,
         estimator="",
         tmle_var=FALSE){
  
  
  cat(paste0(estimator,"\n"))
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("dataset_nums","run_targets_ltmle_simulation2",  "library",
                      "SL.Control","null_sim",
                      "n_bootstrap_samples",
                      "Markov_variables",
                      "gbounds",
                      "n_df",
                      "time","tmle_var"), envir=environment())
  
  clusterEvalQ(cl,lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
    do.call("require", list(X)) 
  }))
  
  
  res=parLapply(cl=cl,dataset_nums, function(z) run_targets_ltmle_simulation2(dataset_num=z,
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

run_targets_ltmle_simulation2 <- function(library="glm",
         SL.Control=NULL,
         dataset_num=NULL,
         n, time=2,
         gbounds=gbounds,
         null_sim=FALSE,
         n_bootstrap_samples=0,
         Markov_variables=NULL,
         tmle_var=tmle_var
){
  
  # if(!is.null(seed)){
  #   set.seed(seed)
  # }
  nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
  nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

  # model <- targets::tar_read_raw("lava_model")
  # simulated_data = simulate_data(lava_model = model, n = n)
  # simulated_data = clean_sim_data(simulated_data, N_time=time)
  simulated_data = readRDS(paste0(here::here(),"/data/sim_data/simulated_data_",dataset_num,".rds"))
  
  if(null_sim){
    simulated_data = data.frame(simulated_data)
    Y = simulated_data[,grepl("dementia",colnames(simulated_data))|grepl("Dead",colnames(simulated_data))|grepl("Censored_",colnames(simulated_data))]
    Y <- Y[sample(nrow(Y)),]
    simulated_data = simulated_data[,!(grepl("dementia",colnames(simulated_data))|grepl("Dead",colnames(simulated_data))|grepl("Censored_",colnames(simulated_data)))]
    simulated_data<- bind_cols(simulated_data, Y)
    simulated_data = as.data.table(simulated_data)
  }
  
  
  # # #set up analysis:
  simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)
  #simulated_data_list$outcome_data
  res = run_ltmle(name_outcome="dementia",
                  time_horizon=time,
                  test=FALSE,
                  outcome_data=simulated_data_list$outcome_data,
                  regimen_data=simulated_data_list$regimen_data,
                  baseline_data=simulated_data_list$sim_baseline_covariates,
                  timevar_data=simulated_data_list$sim_time_covariates,
                  det.Q.function=NULL,# now build-in
                  gbounds=gbounds,
                  B_bootstrap_samples=n_bootstrap_samples,
                  SL.library=library,
                  Markov=Markov_variables,
                  SL.cvControl=SL.Control,
                  tmle_var=tmle_var,
                  verbose=TRUE)
  
  
  res=res[[1]][[1]]
  res=res[1]
  #res
  
  output=summary(res$Ltmle_fit)
  output_iptw= summary(res$Ltmle_fit, estimator="iptw")
  output=bind_rows(output, output_iptw)
  output$dataset_num <- dataset_num
  output
}



ncores=20
estimator="glm"
dataset_nums=c(1,1000)

run_mclapply_targets_ltmle_simulation <- function(ncores=20, estimator="glm", dataset_nums=c(1,1000), library, SL.Control=NULL,gbounds=c(0.01,1)){
  
  reps=dataset_nums[2]/ncores
  full_res=NULL
  
  for(i in dataset_nums[1]:reps){
    start=1+ncores*(i-1)
    stop=i*ncores
    cat(start, "-", stop, "\n")
    res = mclapply_targets_ltmle_simulation2(n_cores=ncores, estimator=estimator,dataset_nums=c(start:stop), library=library, SL.Control=SL.Control, gbounds=gbounds)
    full_res = bind_rows(full_res, res)
  }
  
  return(full_res)  
  
}


# [1] "sim_res_EN"                      "sim_res_Qint_EN"                 "resdf_DetQ_ic_gbound_001"        "resdf_DetQ_ic_gbound_005_9"     
# "sim_res_Qint_ic_glm"              "sim_res_DetQ_ic_v3"             
# [9] "sim_res_DetQ_Qint_ic"            "sim_res_ridge_AUC"               "sim_res_Qint_AUC"                "sim_res_AUC_1se"                
# [13] "sim_res_Qint_AUC_1se"            "sim_res_detQ_rf_v5"              "sim_res_detQ_Qint_rf_v5"         "sim_res_DetQ__ridge_ic_v3"      
# [17] "sim_res_ridge_Qint_v5"    

#sim_res_glm_ic
# res_glm_extra = run_run_mclapply_targets_ltmle_simulation(n_cores=1, estimator="glm",dataset_nums=1, library="glm")

#-------------------------------------------------------------------------------
# Run here
#-------------------------------------------------------------------------------

load(paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))

# system.time({res_glm_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="glm",dataset_nums=c(1,1000), library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))

 # system.time({res_glmnet_extra=run_mclapply_targets_ltmle_simulation(ncores=20, estimator="lasso",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
 # system.time({res_glmnet_undersmooth_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="undersmoothed lasso",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
 # 
 # system.time({res_ridge_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="ridge",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
 # system.time({res_ridge_undersmooth_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="undersmoothed ridge",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
 # system.time({res_EN_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="EN",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
 # system.time({res_EN_undersmooth_extra= run_mclapply_targets_ltmle_simulation(ncores=20, estimator="undersmoothed EN",dataset_nums=c(1,1000), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
 # save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))


# system.time({res_glm_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=20, estimator="glm",dataset_nums=c(1,1000), gbounds=c(0,1), library="glm")})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_glmnet_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=20, estimator="lasso",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_glmnet_undersmooth_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=20, estimator="undersmoothed lasso",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=1))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_ridge_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=50, estimator="ridge",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_ridge_undersmooth_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=50, estimator="undersmoothed ridge",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_EN_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=50, estimator="EN",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="min_lambda",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# system.time({res_EN_undersmooth_untruncated_extra=run_mclapply_targets_ltmle_simulation(ncores=50, estimator="undersmoothed EN",dataset_nums=c(1,1000), gbounds=c(0,1), library="glmnet", SL.Control=list(selector="undersmooth",alpha=0.5))})
# save(list=ls(pattern = "res_"), file=paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))



