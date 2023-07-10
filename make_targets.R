
library(here)
library(targets)
source(paste0(here::here(),"/_targets.R"))
tar_make(garbage_collection = TRUE)



res1 = tar_read(test_results_1)
res2 = tar_read(test_results_2)
res3 = tar_read(test_results_3)
res_boot1 = tar_read(test_results_boot_2)

summary(res1)
summary(res2)
summary(res3)
res_boot1


res_boot = tar_read(glmnet_res_boot)
saveRDS(res_boot, file = paste0(here::here(),"/scratch/res_boot_v1.RDS"))

res_boot <- res_boot %>% filter(Target_parameter=="ATE", Estimator=="tmle")
boot_CIs <- res_boot %>% group_by(iteration) %>%
  summarise(
    mean_est = mean(estimate),
    CI1=quantile(estimate,.025),
    CI2=quantile(estimate,.975)
  )
mean(boot_CIs$mean_est)
CI_widths = boot_CIs$CI2 - boot_CIs$CI1
mean(CI_widths)
median(CI_widths)

res = tar_read(glmnet_res)
saveRDS(res, file = paste0(here::here(),"/scratch/res_v1.RDS"))
res <- res %>% filter(Target_parameter=="ATE", Estimator=="tmle")
mean(res$estimate)
mean(res$upper - res$lower)


#bootstrap is slightly lower with updated res


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
# save the current results
# update bootstrap to use L2 to pred Y2 instead of Y3, etc.
#calculate current coverage
#test IC with and without Markov process
# Simulate the really simple data using code from other repo to check
# do simulation without death or longitudinal confounding
# do simulation without death
#hand calc true RR just using coefficients
#update truth function to calculate average across many seeds

# event_dementia_3 using 94989 observations.
# Estimate: framing formula Q.kplus1 ~ baseline_variables + insulin_2 + any.malignancy_2 + chronic.pulmonary.disease_2 + hypertension_2 + myocardial.infarction_2 + ischemic.heart.disease_2 + heart.failure_2 + renal.disease_2 + sglt2_inhib_2 + glp1_0 + glp1_1 + glp1_2 into Y and X...


df=tar_read(sim_data)[[1]]
fit=run_Ltmle(d=df,
                   time_horizon=4,
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
              concurrentY=TRUE,
              SL.cvControl=list(selector="min_lambda",alpha=0),
              #SL.library = "glm",
                   #deterministic.Q.function = det.Q.function,
                   #test = FALSE,
                  Markov_vars=NULL)

#NOTE! Need to incorporate A3 correctly in predicting Y3, etc
#need to fix censoring formula


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

