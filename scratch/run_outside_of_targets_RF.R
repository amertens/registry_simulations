
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
seeds1=c(136632L, 702347L, 31157L, 652859L, 501151L, 344524L, 672750L, 
         247494L, 405801L, 994464L, 590739L, 952491L, 20091L, 310026L, 
         193136L, 401333L, 291832L, 525393L, 614165L, 520503L, 574780L, 
         593146L, 964646L, 240460L, 80989L, 750527L, 520317L, 875104L, 
         147399L, 81091L, 320258L, 685634L, 377406L, 765353L, 340573L, 
         56025L, 737603L, 284593L, 291591L, 184969L, 712818L, 311320L, 
         151148L, 361101L, 796200L, 647711L, 2788L, 415658L, 701673L, 
         461213L)
set.seed(12345)
seeds2=sample(1000001:2000000, 250, replace=FALSE)
set.seed(12345)
seeds3=sample(2000001:3000000, 200, replace=FALSE)
seeds_rf=c(seeds1, seeds2, seeds3)

#-------------------------------------------------------
# Custom RF function
#-------------------------------------------------------

SL.ranger.custom<-function (Y, X, newX, family, obsWeights, num.trees = 100, mtry = floor(sqrt(ncol(X))), 
                            write.forest = TRUE, probability = family$family == "binomial", 
                            min.node.size = ifelse(family$family == "gaussian", 
                                                   5, 1), replace = TRUE, sample.fraction = ifelse(replace, 
                                                                                                   1, 0.632), 
                            #number of cores to use:
                            num.threads = 10, 
                            verbose = T, ...) 
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

rf_res_list <- vector("list", length(seeds_rf))


#-------------------------------------------------------
# Custom RF function
#-------------------------------------------------------


for(i in 1:seeds_rf){
    
    res=NULL
    system.time({res=run_targets_ltmle_simulation(seed=seeds_rf[i], library="SL.ranger.custom",
                                                                        SL.Control=NULL,
                                                                        n_bootstrap_samples=0,
                                                                        Markov_variables=Markov_variables,
                                                                        gbounds=c(0.01,1),
                                                                        n=100000,
                                                                        time=10)})
    res$estimator="random forest markov"
    rf_res_list[[i]] <- res
    saveRDS(rf_res_list,  file=paste0(here::here(),"/data/sim_results/sim_res_RF.RDS"))
}



