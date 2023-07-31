


# library="glmnet"
# SL.Control=list(selector="undersmooth",alpha=1)
# n=n_df
# time=10
# n_bootstrap_samples=0
# rep_seed=12345

run_targets_ltmle_simulation_comparison <- function(library="glm",
                                                    rep_seed=12345,
                                               SL.Control=NULL,
                                               n, 
                                               time=2,
                                               n_bootstrap_samples=0,
                                               Markov_variables=NULL){
  
  set.seed(rep_seed)
  
  # source("data/coefs.txt")
  # model= get_lava_model(time_horizon = time, coefs = coefs)
  model <- targets::tar_read_raw("lava_model")
  simulated_data = simulate_data(lava_model = model, n = n)
  simulated_data = clean_sim_data(simulated_data, N_time=time)
     
    # # #set up analysis:
     simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)

    res1 = run_ltmle(name_outcome="dementia",
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
    
    
    res1=res1[[1]][[1]]
    res1=res1[1]
    res1
    

}

