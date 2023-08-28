
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(ranger)
library(parallel)

lapply(c("ranger","fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})


# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)


# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
#set
#Longitudinal variables possibly following the markov process
Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                   "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                   "bb","ccb","rasi","thiazid",
                   "loop","mra","copd_med"  )

#seedlists
set.seed(12345)
seeds_null=sample(0:1000000, 500, replace=FALSE)

#-------------------------------------------------------
# Custom RF function
#-------------------------------------------------------

SL.ranger.custom<-function (Y, X, newX, family, obsWeights, num.trees = 10, mtry = floor(sqrt(ncol(X))), 
                            write.forest = TRUE, probability = family$family == "binomial", 
                            min.node.size = ifelse(family$family == "gaussian", 
                                                   5, 1), replace = TRUE, sample.fraction = ifelse(replace, 
                                                                                                   1, 0.632), 
                            #number of cores to use:
                            num.threads = 90, 
                            verbose = F, ...) 
{
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, 
                                               X), num.trees = num.trees, mtry = mtry, min.node.size = min.node.size, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        probability = probability, num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}

rf_res_list <- vector("list", length(seeds_null))
rf_res_list <- readRDS(paste0(here::here(),"/data/sim_results/sim_res_RF_NULL.RDS"))

#-------------------------------------------------------
# Custom RF function
#-------------------------------------------------------


for(i in 1:seeds_null){
    
  cat(i,"\n")
    res=NULL
    system.time({res=run_targets_ltmle_simulation(null_sim=TRUE, seed=seeds_null[i], library="SL.ranger.custom",
                                                                        SL.Control=NULL,
                                                                        n_bootstrap_samples=0,
                                                                        Markov_variables=Markov_variables,
                                                                        gbounds=c(0.01,1),
                                                                        n=100000,
                                                                        time=10)})
    res$estimator="random forest markov"
    rf_res_list[[i]] <- res
    saveRDS(rf_res_list,  file=paste0(here::here(),"/data/sim_results/sim_res_RF_NULL.RDS"))
}





