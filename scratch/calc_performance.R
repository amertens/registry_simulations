
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

load(paste0(here::here(),"/data/sim_results/sim_res_seeds1.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))

res1=mget(ls(pattern = "res_"))

res1= data.table::rbindlist(l=res1, use.names=TRUE, fill=TRUE, idcol="analysis")
res1$analysis<-gsub("_[0-9]+$","",res1$analysis)
table(res1$analysis)
# 
# res2=list( tar_read(res_glm_1),
#              tar_read(res_glmnet_1),
#              tar_read(res_glmnet_undersmooth_1),
#              tar_read(res_ridge_1),
#              tar_read(res_ridge_undersmooth_1),
#              tar_read(res_EN_1),
#              tar_read(res_EN_undersmooth_1),
#              tar_read(res_glm_markov_1),
#              tar_read(res_glmnet_markov_1),
#              #tar_read(res_glmnet_undersmooth_markov_1),
#              #tar_read(res_ridge_markov_1),
#              #tar_read(res_ridge_undersmooth_markov_1),
#              tar_read(res_EN_markov_1),
#              tar_read(res_EN_undersmooth_markov_1),
#              tar_read(res_glm_untruncated_1)#,
#              #tar_read(res_glmnet_untruncated_1),
#              #tar_read(res_glmnet_undersmooth_untruncated_1),
#              # tar_read(res_ridge_untruncated_1),
#              # tar_read(res_ridge_undersmooth_untruncated_1),
#              # tar_read(res_EN_untruncated_1),
#              # tar_read(res_EN_undersmooth_untruncated_1),
#              # tar_read(res_glm_markov_untruncated_1),
#              # tar_read(res_glmnet_markov_untruncated_1),
#              # tar_read(res_glmnet_undersmooth_markov_untruncated_1),
#              # tar_read(res_ridge_markov_untruncated_1),
#              # tar_read(res_ridge_undersmooth_markov_untruncated_1),
#              # tar_read(res_EN_markov_untruncated_1),
#              # tar_read(res_EN_undersmooth_markov_untruncated_1)
#              )
# names(res2) <- c("res_glm_tr",
#                  "res_glmnet_tr",
#                  "res_glmnet_undersmooth_tr",
#                  "res_ridge_tr",
#                  "res_ridge_undersmooth_tr",
#                  "res_EN_tr",
#                  "res_EN_undersmooth_tr",
#                  "res_glm_markov_tr",
#                  "res_glmnet_markov_tr",
#                  "res_EN_markov_tr",
#                  "res_EN_undersmooth_markov_tr",
#                  "res_glm_untruncated_tr")
# res2= data.table::rbindlist(l=res2, use.names=TRUE, fill=TRUE, idcol="analysis")
# table(res2$analysis)
# 
# 
# res=bind_rows(res1, res2)
res=res1
#res=c(res1,res2)
# drop <- rep(FALSE, length(res))
# for(i in 1:length(res)){
#   if(class(res[[i]])[1]!="data.table"){
#     drop[i] <- TRUE
#   }
# }
# 
# res[drop] <- NULL
truth=tar_read(truth)
time=10

res$estimator = gsub("res_","",res$analysis)
res$estimator = gsub("_tr","",res$estimator)
res$estimator = gsub("_[0-9]+$","",res$estimator)
table(res$estimator)
sim_perf_tab = calc_sim_performance(
  res=res,
  truth=tar_read(truth),
  time=10)
sim_perf_tab




