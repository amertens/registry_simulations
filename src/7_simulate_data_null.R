
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(data.table)
library(tidyverse)
library(parallel)
library(lava)


# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)


source("data/coefs.txt")

coefs_null=coefs
coefs_null$outcome_coef$dementia_1[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_1))] <- 0
coefs_null$outcome_coef$dementia_2[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_2))] <- 0
coefs_null$outcome_coef$dementia_3[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_3))] <- 0
coefs_null$outcome_coef$dementia_4[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_4))] <- 0
coefs_null$outcome_coef$dementia_5[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_5))] <- 0
coefs_null$outcome_coef$dementia_6[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_6))] <- 0
coefs_null$outcome_coef$dementia_7[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_7))] <- 0
coefs_null$outcome_coef$dementia_8[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_8))] <- 0
coefs_null$outcome_coef$dementia_9[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_9))] <- 0
coefs_null$outcome_coef$dementia_10[grepl("GLP1RA",names(coefs_null$outcome_coef$dementia_10))] <- 0


coefs_null$comp.event_coef$Dead_1[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_1))] <- 0
coefs_null$comp.event_coef$Dead_2[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_2))] <- 0
coefs_null$comp.event_coef$Dead_3[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_3))] <- 0
coefs_null$comp.event_coef$Dead_4[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_4))] <- 0
coefs_null$comp.event_coef$Dead_5[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_5))] <- 0
coefs_null$comp.event_coef$Dead_6[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_6))] <- 0
coefs_null$comp.event_coef$Dead_7[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_7))] <- 0
coefs_null$comp.event_coef$Dead_8[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_8))] <- 0
coefs_null$comp.event_coef$Dead_9[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_9))] <- 0
coefs_null$comp.event_coef$Dead_10[grepl("GLP1RA",names(coefs_null$comp.event_coef$Dead_10))] <- 0

model_null = get_lava_model(time_horizon = 10, coefs = coefs_null)


#Check null simulation by checking the truth
set.seed(123)
modelA1 = get_lava_model(time_horizon = 10, coefs = coefs_null, counterfactual_A=1)
modelA0 = get_lava_model(time_horizon = 10, coefs = coefs_null, counterfactual_A=0)

d.always = simulate_data(lava_model = modelA1, n = 1000000)
d.never = simulate_data(lava_model = modelA0, n = 1000000)

Ya0_t10 <- mean(d.never$dementia_10,na.rm=T)
Ya1_t10 <- mean(d.always$dementia_10,na.rm=T)
Ya1_t10 - Ya0_t10
Ya1_t10 / Ya0_t10

Ya0_t10 <- mean(d.never$Dead_9,na.rm=T)
Ya1_t10 <- mean(d.always$Dead_9,na.rm=T)
Ya1_t10 - Ya0_t10
Ya1_t10 / Ya0_t10


# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
#Longitudinal variables possibly following the markov process
Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                   "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                   "bb","ccb","rasi","thiazid",
                   "loop","mra","copd_med"  )

#set seed
seed=3456
simulated_data_list =NULL
i=1
for(i in 367:1000){
  cat(i, "\n")
  set.seed(seed+i)
  simulated_data = simulate_data(lava_model = model_null, n = n_df)
  saveRDS(simulated_data, file = paste0(here::here(),"/data/sim_data/null/simulated_data_",i,".rds"))
}

