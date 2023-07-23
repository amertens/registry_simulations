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

run_targets_ltmle_simulation_bootstrap <- function(library="glmnet",
                                               SL.Control=NULL,
                                               n, time=2,
                                               B_bootstrap_samples=0,
                                               Markov_variables=NULL){
  
  browser()

  model <- targets::tar_read_raw("lava_model")
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data <- data.table::as.data.table(simulated_data)
  dboot <- simulated_data[sample(.N, nrow(simulated_data),replace=TRUE)]
  #hack: remake unique PNR's for merge
  dboot$pnr <- 1:nrow(dboot)
  dboot
  simulated_data = clean_sim_data(dboot, N_time=time)

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
                                             SL.library=library,
                                             Markov=Markov_variables,
                                             SL.cvControl=SL.Control,
                                             verbose=TRUE)


    res=res[[1]][[1]]
    res=res[1]
    res



}

