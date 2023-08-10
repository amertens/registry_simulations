### run_ltmle_simulation.R --- 
#----------------------------------------------------------------------
## Author: Andrew Mertens
## Created: Jul 17 2023 (14:14) 
## Version: 
## Last-Updated: Jul 17 2023 (15:23) 
##           By: Andrew Mertens
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log: 
#----------------------------------------------------------------------
## 
### Code:

#Wrapper functions for simulating data then analyzing it 
#(ideally across many estimation hyperparameters)

# seed=1234
# library="glmnet"
# SL.Control=list(selector="undersmooth",alpha=1)
# n=10000
# time=2
# Markov_variables=Markov_variables
# n_bootstrap_samples=0




run_targets_ltmle_simulation <- function(library="glm",
                                               SL.Control=NULL,
                                               seed=NULL,
                                               n, time=2,
                                               gbounds=gbounds,
                                               null_sim=FALSE,
                                               n_bootstrap_samples=0,
                                               Markov_variables=NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
  nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
  
  # source("data/coefs.txt")
  # model= get_lava_model(time_horizon = time, coefs = coefs)
  model <- targets::tar_read_raw("lava_model")
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data = clean_sim_data(simulated_data, N_time=time)
  
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
                                             verbose=TRUE)
    
    
    res=res[[1]][[1]]
    res=res[1]
    #res

    output=summary(res$Ltmle_fit)
    output
}

