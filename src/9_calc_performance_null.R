
rm(list=ls())
library(tidyverse)
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

load(paste0(here::here(),"/data/sim_results/sim_res_null_markov_update_v2.Rdata"))
# res_RF_null <- readRDS(paste0(here::here(),"/data/sim_results/sim_res_RF_NULL.RDS"))
# "C:\Users\andre\Documents\jici\registry_simulations\data\sim_results\sim_res_RF_NULL.RDS"

res_EN2=res_EN                            
res_EN_undersmooth2=res_EN_undersmooth                
res_EN_undersmooth_untruncated2=res_EN_undersmooth_untruncated     
res_EN_untruncated2=res_EN_untruncated                
res_glm2=res_glm                           
res_glm_untruncated2=res_glm_untruncated               
res_glmnet2=res_glmnet                         
res_glmnet_undersmooth2=res_glmnet_undersmooth            
res_glmnet_undersmooth_untruncated2=res_glmnet_undersmooth_untruncated 
res_glmnet_untruncated2=res_glmnet_untruncated            
res_ridge2=res_ridge                          
res_ridge_undersmooth2=res_ridge_undersmooth             
res_ridge_undersmooth2_untruncated=res_ridge_undersmooth_untruncated  
res_ridge_untruncated2=res_ridge_untruncated   

load(paste0(here::here(),"/data/sim_results/sim_res_null_markov_update_v3.Rdata"))


res=mget(ls(pattern = "res_"))

res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")
res$estimator<-gsub("_[0-9]+$","",res$analysis)
res$estimator = gsub("res_","",res$estimator)
res$estimator = gsub("_tr","",res$estimator)
res$estimator = gsub("2","",res$estimator)
res$analysis  = gsub("2","",res$analysis )
table(res$estimator)

saveRDS(res, file=paste0(here::here(),"/data/sim_results/sim_res_null.rds"))

truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds")) %>% subset(., select=-c(meanYa0))
truth[,3]=truth[,6]=1
truth[,4]=truth[,7]=0

sim_perf_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)
sim_perf_tab

write.csv(sim_perf_tab, paste0(here::here(),"/data/sim_perf_null.csv"))
