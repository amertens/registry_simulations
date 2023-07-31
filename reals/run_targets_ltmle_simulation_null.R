### run_ltmle_simulation_null.R --- 
#----------------------------------------------------------------------
## Author: Andrew Mertens
## Created: Jul 26 2023 (14:14) 
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

run_targets_ltmle_simulation_null <- function(library="glm",
                                               SL.Control=NULL,
                                               n=100000, 
                                               time=2,
                                               n_bootstrap_samples=0,
                                               Markov_variables=NULL){
  
  
  model <- targets::tar_read_raw("lava_model")
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data = clean_sim_data(simulated_data, N_time=time)
  simulated_data = data.frame(simulated_data)
  Y = simulated_data[,grepl("dementia",colnames(simulated_data))|grepl("Dead",colnames(simulated_data))|grepl("Censored_",colnames(simulated_data))]
  Y <- Y[sample(nrow(Y)),]
  simulated_data = simulated_data[,!(grepl("dementia",colnames(simulated_data))|grepl("Dead",colnames(simulated_data))|grepl("Censored_",colnames(simulated_data)))]
  simulated_data<- bind_cols(simulated_data, Y)
  simulated_data = as.data.table(simulated_data)

  # #set up analysis:
  simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)

  res = run_ltmle(name_outcome="dementia",
                  time_horizon=time,
                  test=FALSE,
                  outcome_data=simulated_data_list$outcome_data,
                  regimen_data=simulated_data_list$regimen_data,
                  baseline_data=simulated_data_list$sim_baseline_covariates,
                  timevar_data=simulated_data_list$sim_time_covariates,
                  det.Q.function=NULL,# now build-in
                  B_bootstrap_samples=n_bootstrap_samples,
                  SL.library=library,
                  Markov=Markov_variables,
                  SL.cvControl=SL.Control,
                  verbose=TRUE)


  res=res[[1]][[1]]
  res=res[1]
  res
  
}

