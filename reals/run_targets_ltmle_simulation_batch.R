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

run_targets_ltmle_simulation_batch <- function(library="glm",
                                               SL.Control=NULL,
                                               n, time=2,
                                               Markov_variables=NULL){
  
  

  
  source("data/coefs.txt")
  model= get_lava_model(time_horizon = time, coefs = coefs)
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data = clean_sim_data(simulated_data, N_time=time)
     
    # # #set up analysis:
     simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)
     #simulated_data_list$outcome_data
    estimate = run_ltmle(name_outcome="dementia",
                                             time_horizon=time,
                                             test=FALSE,
                                             outcome_data=simulated_data_list$outcome_data,
                                             regimen_data=simulated_data_list$regimen_data,
                                             baseline_data=simulated_data_list$sim_baseline_covariates,
                                             timevar_data=simulated_data_list$sim_time_covariates,
                                             det.Q.function=NULL,# now build-in
                                             SL.library=library,
                                             Markov=Markov_variables,
                                             SL.cvControl=SL.Control,
                                             verbose=TRUE)
    return(estimate)

}

