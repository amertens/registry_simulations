
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
load(paste0(here::here(),"/data/sim_results/sim_res_markov_alt (3).Rdata"))


ls(pattern = "res_")


res=mget(ls(pattern = "res_"))

res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")


res
res %>% filter(Target_parameter=="ATE")  %>% group_by(analysis, Estimator) %>% mutate(o.var=var(estimate)) %>%
  summarize(abs_bias= mean(abs(estimate - (-0.0071862 ))), 
            o.coverage=mean(estimate-1.96*sd(estimate)< -0.0071862 & -0.0071862 < estimate+1.96*sd(estimate))*100,
            coverage= mean(  lower < -0.0071862 &     upper > -0.0071862 ))


table(res$estimator)
table(res$analysis)
#TEMP
res$estimator <- res$analysis


saveRDS(res, file=paste0(here::here(),"/data/sim_results/sim_res.rds"))


truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds")) %>% subset(., select=-c(meanYa0))
#truth_old =tar_read(truth)

sim_perf_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10, 
  mean=FALSE)
sim_perf_tab %>% filter(Estimator=="tmle") %>% arrange(O_coverage_RD)
sim_perf_tab %>% filter(Estimator=="tmle") %>% arrange(O_coverage_RD)

sim_perf_tab[sim_perf_tab$estimator=="res_ridge_markov_undersmooth" & sim_perf_tab$Estimator=="tmle",]

write.csv(sim_perf_tab, paste0(here::here(),"/data/sim_perf_1000reps.csv"))

#look for smallest variance with 1% of 95% oracle coverage
sim_perf_tab %>% filter(abs(O_coverage_RD-95) <= 1) %>% arrange(estimator_variance_RD  )


trueRD= -0.0071862
trueRD= -0.007200 #median

trueRD=-0.0071820  
trueRD=-0.007195

trueRD=-0.0071970
trueRD=-0.007220
res %>% filter(Target_parameter=="ATE") %>%
  group_by(Estimator, estimator) %>% 
  summarise(abs_bias_RD=mean(abs(estimate-trueRD)),
            #estimator_variance_RD=mean(((estimate)-mean((estimate)))^2),
            mean_variance_RD=mean((std.err)^2),
            #bias_se_ratio_RD=abs_bias_RD/sqrt(mean_variance_RD),
            #coverage_RD=mean(lower<=trueRD & trueRD<=upper)*100,
            O_coverage_RD=mean(estimate-1.96*sd(estimate)< trueRD & trueRD < estimate+1.96*sd(estimate))*100) %>%
  filter(estimator %in% c("res_glm","res_glm_untruncated","res_ridge_undersmooth","res_ridge_undersmooth_untruncated"), Estimator!="iptw") %>%
  arrange(mean_variance_RD)
