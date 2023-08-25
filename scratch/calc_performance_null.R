
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

list.files(paste0(here::here(),"/data/sim_results/"))
# load(paste0(here::here(),"/data/sim_results/sim_res_full.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seeds_null.Rdata"))


res=mget(ls(pattern = "res_"))

res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")
res$estimator<-gsub("_[0-9]+$","",res$analysis)
res$estimator = gsub("res_","",res$estimator)
res$estimator = gsub("_tr","",res$estimator)

res <- res %>% group_by(estimator) %>% slice(1:4000)
table(res$estimator)

saveRDS(res, file=paste0(here::here(),"/data/sim_results/sim_res_null.rds"))

#TO DO
#add boxplots, pick estimator, run bootstrap
#add the IPTW estimates

# truth=tar_read(truth)
# saveRDS(truth, file=paste0(here::here(),"/data/sim_results/truth.rds"))
truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds"))

truth=tar_read(truth)
truth[10,3]=1
truth[10,4]=0

sim_perf_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)
sim_perf_tab

write.csv(sim_perf_tab, paste0(here::here(),"/data/sim_perf_null.csv"))
