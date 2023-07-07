
library(here)
library(targets)
source(paste0(here::here(),"/_targets.R"))
tar_make(garbage_collection = TRUE)

res_boot = tar_read(test_results_bootstrap)


res1 = tar_read(test_results_1)
res2 = tar_read(glm_res)
res = res2 %>% filter(Target_parameter=="ATE", Estimator=="tmle")
mean(res$estimate)
res = res2 %>% filter(Target_parameter=="RelativeRisk", Estimator=="tmle")
exp(mean(log(res$estimate)))


res3 = tar_read(glmnet_res)
res = res3 %>% filter(Target_parameter=="ATE", Estimator=="tmle")
mean(res$estimate)
res = res3 %>% filter(Target_parameter=="RelativeRisk", Estimator=="tmle")
exp(mean(log(res$estimate)))

truth = tar_read(truth)
truth$truth_df

#Scratch
#To do next: 
#make a bootstrapping option
# Simulate the really simple data using code from other repo to check
# do simulation without death or longitudinal confounding
# do simulation without death
#hand calc true RR just using coefficients
#update truth function to calculate average across many seeds

# event_dementia_3 using 94989 observations.
# Estimate: framing formula Q.kplus1 ~ baseline_variables + insulin_2 + any.malignancy_2 + chronic.pulmonary.disease_2 + hypertension_2 + myocardial.infarction_2 + ischemic.heart.disease_2 + heart.failure_2 + renal.disease_2 + sglt2_inhib_2 + glp1_0 + glp1_1 + glp1_2 into Y and X...


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
test=run_ltmle_sim(sim_d_list=sim_data, time_horizon=11)
head(test)
table(test$iteration)


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

