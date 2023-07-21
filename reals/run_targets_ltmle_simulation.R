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

run_targets_ltmle_simulation <- function(lava_model = lava_model, n=10000, library="glm", time=2){


    simulated_data <- simulate_data(lava_model = lava_model, n = n)
    #set up analysis:
    simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)
    estimate = run_ltmle(name_outcome="dementia",
                                             time_horizon=time,
                                             test=FALSE,
                                             outcome_data=simulated_data_list$outcome_data,
                                             regimen_data=simulated_data_list$regimen_data,
                                             baseline_data=simulated_data_list$sim_baseline_covariates,
                                             timevar_data=simulated_data_list$sim_time_covariates,
                                             det.Q.function=NULL,# now build-in
                                             SL.library=library,
                                             #SL.cvControl=list(selector="undersmooth",alpha=1),
                                             verbose=TRUE)
    res=summary(estimate[[1]][[1]]$Ltmle_fit)
    
    res.iptw = summary(estimate[[1]][[1]]$Ltmle_fit,"iptw")
    res=bind_rows(res, res.iptw)
    # #add in
    # # res$Qform = list(fit$formulas$Qform)
    # # res$gform = list(fit$formulas$gform)
    # res$analysis_name=analysis_name
    # res$iteration=iteration
    res=data.frame(res)
    return(res)
}

