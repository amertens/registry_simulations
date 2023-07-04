
library(here)
library(targets)
source(paste0(here::here(),"/_targets.R"))
tar_make()


#Scratch
#To do next: 
#make a bootstrapping option
# Simulate the really simple data using code from other repo to check
# do simulation without death or longitudinal confounding
# do simulation without death
#hand calc true RR just using coefficients
#update truth function to calculate average across many seeds

df=tar_read(sim_data)[[1]]
fit=run_Ltmle(d=df,
                   #time_horizon=3,
                   # name_outcome = "event_dementia",
                   # name_regimen = "glp1",
                   # name_censoring = "censor",
                   # censored_label = "censored",
                   # name_comp.event = "event_death",
                   #baseline_vars=baseline_vars,
                   #long_covariates=long_covariates,
                  # treatment_vars="glp1",
                  # outcome_vars=c("event_dementia","censor","event_death"),
              SL.library = "glmnet",
              SL.cvControl=list(selector="min_lambda",alpha=0),
              #SL.library = "glm",
                   #deterministic.Q.function = det.Q.function,
                   #test = FALSE,
                  Markov_vars=Markov_variables)
summary(fit)


sim_data=tar_read(sim_data)
test=run_ltmle_sim(sim_d_list=sim_data, time_horizon=2)
head(test)



#----------
tar_visnetwork(targets_only = TRUE)


targets::tar_meta(fields = warnings, complete_only = TRUE)


#To do:
#either create or load the data with each iteration of the simulation, rather than loading a list into data
#set up parallelization within the targets system
#figure out how to use Thomas's glmnet updates

#Compare bootstrap and IC CI's (first with glm, then glmnet)
#Get targets based sim results for simple sim (short followup, no death)

#update clean_ltmle_res function to format the results


#check
d <- tar_read(sim_data)
d[[1]]

d <- tar_read(test_results_2)
d

