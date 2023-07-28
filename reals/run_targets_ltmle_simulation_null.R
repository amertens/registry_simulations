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
                                               n, time=2,
                                               n_bootstrap_samples=0,
                                               Markov_variables=NULL){
  
  # source("data/coefs.txt")
  # model= get_lava_model(time_horizon = time, coefs = coefs)
  model <- targets::tar_read_raw("lava_model")
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data = clean_sim_data(simulated_data, N_time=time)
    Y <- simulated_data %>% select(starts_with("dementia")|starts_with("dead")|starts_with("censor_"))
    Y <- Y[sample(nrow(Y)),]
    simulated_data <- simulated_data %>% select(!starts_with("dementia")&!starts_with("dead")&!starts_with("censor_"))
    dsimulated_data<- bind_cols(simulated_data, Y)
  
  
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
                  B_bootstrap_samples=n_bootstrap_samples,
                  SL.library=library,
                  Markov=Markov_variables,
                  SL.cvControl=SL.Control,
                  verbose=TRUE)
  
  
  res=res[[1]][[1]]
  res=res[1]
  res
  
  # res=summary(res[grepl("Ltmle_fit",names(res))])
  # res
  #summary(res[[1]])
  #(res[grepl("Ltmle_fit",names(res))])
  #names(res[[1]][[1]])
  # Ltmle_fit=NULL
  # try(Ltmle_fit = summary(res[[1]][[1]]$Ltmle_fit))
  # Ltmle_fit
  #return(res_tab)
  
}

