
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
load(paste0(here::here(),"/data/sim_results/sim_res_markov_extra.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_markov_RF.Rdata"))

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

load(paste0(here::here(),"/data/sim_results/sim_res_markov_update.Rdata"))

res_EN=bind_rows(res_EN2, res_EN) %>% slice(1:8000)
res_EN_undersmooth=bind_rows(res_EN_undersmooth2, res_EN_undersmooth) %>% slice(1:8000)
res_EN_undersmooth_untruncated=bind_rows(res_EN_undersmooth_untruncated2, res_EN_undersmooth_untruncated) %>% slice(1:8000)
res_EN_untruncated=bind_rows(res_EN_untruncated2, res_EN_untruncated) %>% slice(1:8000)
res_glm=bind_rows(res_glm2, res_glm) %>% slice(1:8000)
res_glm_untruncated=bind_rows(res_glm_untruncated2[1:8,], res_glm_untruncated) %>% slice(1:8000)
res_glmnet=bind_rows(res_glmnet2, res_glmnet) %>% slice(1:8000)
res_glmnet_undersmooth=bind_rows(res_glmnet_undersmooth2, res_glmnet_undersmooth) %>% slice(1:8000)
res_glmnet_undersmooth_untruncated=bind_rows(res_glmnet_undersmooth_untruncated2, res_glmnet_undersmooth_untruncated) %>% slice(1:8000)
res_glmnet_untruncated=bind_rows(res_glmnet_untruncated2, res_glmnet_untruncated) %>% slice(1:8000)
res_ridge=bind_rows(res_ridge2, res_ridge) %>% slice(1:8000)
res_ridge_undersmooth=bind_rows(res_ridge_undersmooth2, res_ridge_undersmooth) %>% slice(1:8000)
res_ridge_undersmooth_untruncated=bind_rows(res_ridge_undersmooth2_untruncated, res_ridge_undersmooth_untruncated) %>% slice(1:8000)
res_ridge_untruncated=bind_rows(res_ridge_untruncated2, res_ridge_untruncated) %>% slice(1:8000)

rm(res_EN2, res_EN_undersmooth2, res_EN_undersmooth_untruncated2, res_EN_untruncated2, res_glm2, res_glm_untruncated2, res_glmnet2, res_glmnet_undersmooth2, res_glmnet_undersmooth_untruncated2, res_glmnet_untruncated2, res_ridge2, res_ridge_undersmooth2, res_ridge_undersmooth2_untruncated, res_ridge_untruncated2)


ls(pattern = "res_")

which(!(1:1000 %in% res_ridge_undersmooth_untruncated$dataset_num))
which(!(1:1000 %in% res_glmnet_undersmooth_untruncated$dataset_num))


res_ridge_undersmooth_untruncated %>% group_by(Estimator) %>% 
  summarize(N=n(), mean(estimate), mean(std.err))

res_ridge_undersmooth %>% filter(Target_parameter=="ATE")  %>% group_by(Estimator) %>% mutate(o.var=var(estimate)) %>%
  summarize(mean_est=mean(estimate),
            abs_bias= mean(abs(estimate - (-0.0071862 ))),
            sd(estimate),
            o.coverage.lb=mean(estimate-1.96*sd(estimate)),
            o.coverage.ub=mean(estimate+1.96*sd(estimate)),
            o.coverage=mean(estimate-1.96*sd(estimate)< -0.0071862 & -0.0071862 < estimate+1.96*sd(estimate))*100,
            coverage= mean(  lower < -0.0071862 &     upper > -0.0071862 ))

res_glm %>% filter(Target_parameter=="ATE")  %>% group_by(Estimator) %>% mutate(o.var=var(estimate)) %>%
  summarize(mean_est=mean(estimate),
            abs_bias= mean(abs(estimate - (-0.0071862 ))),
            sd(estimate),
            o.coverage.lb=mean(estimate-1.96*sd(estimate)),
            o.coverage.ub=mean(estimate+1.96*sd(estimate)),
            o.coverage=mean(estimate-1.96*sd(estimate)< -0.0071862 & -0.0071862 < estimate+1.96*sd(estimate))*100,
            coverage= mean(  lower < -0.0071862 &     upper > -0.0071862 ))

res_ridge_undersmooth_untruncated %>% filter(Target_parameter=="ATE")  %>% group_by(Estimator) %>% mutate(o.var=var(estimate)) %>%
  summarize(abs_bias= mean(abs(estimate - (-0.0071862 ))), 
            o.coverage=mean(estimate-1.96*sd(estimate)< -0.0071862 & -0.0071862 < estimate+1.96*sd(estimate))*100,
            coverage= mean(  lower < -0.0071862 &     upper > -0.0071862 ))

res_glmnet_undersmooth_untruncated %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))
res_glm %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))
res_RF %>% group_by(Estimator) %>% 
  summarize( mean(estimate), mean(std.err))


rm(res_glmnet2 )
rm(res_test )
res=mget(ls(pattern = "res_"))

res= data.table::rbindlist(l=res, use.names=TRUE, fill=TRUE, idcol="analysis")

res %>% filter(Target_parameter=="ATE")  %>% group_by(Estimator) %>% mutate(o.var=var(estimate)) %>%
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
