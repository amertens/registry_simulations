
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)

lapply(c("targets","fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

list.files(paste0(here::here(),"/data/sim_results/"))
# load(paste0(here::here(),"/data/sim_results/sim_res_extra.Rdata"))
# load(paste0(here::here(),"/data/sim_results/sim_res_full.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_markov.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_markov_RF.Rdata"))

ls(pattern = "res_")

res_EN_markov %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))


res_glmnet_markov_undersmooth_untruncated %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))
res_glm_markov %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))
res_random_forest %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))



res=mget(ls(pattern = "res_"))

res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")

table(res$estimator)
table(res$analysis)
#TEMP
res$estimator <- res$analysis


saveRDS(res, file=paste0(here::here(),"/data/sim_results/sim_res.rds"))


truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds"))
#truth_old =tar_read(truth)

sim_perf_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)
sim_perf_tab %>% filter(Estimator=="tmle") %>% arrange(O_coverage_RD)

sim_perf_tab[sim_perf_tab$estimator=="res_ridge_markov_undersmooth" & sim_perf_tab$Estimator=="tmle",]

write.csv(sim_perf_tab, paste0(here::here(),"/data/sim_perf_1000reps.csv"))

#look for smallest variance with 1% of 95% oracle coverage
sim_perf_tab %>% filter(abs(O_coverage_RD-95) <= 1) %>% arrange(mean_variance_RD )

