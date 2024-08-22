

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
load(paste0(here::here(),"/data/sim_results/sim_res_seeds1.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed2.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed3.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed4.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_seed5.Rdata"))
load(paste0(here::here(),"/data/sim_results/sim_res_markov.Rdata"))


res <- bind_rows(res_glmnet_markov_undersmooth, res_glmnet_undersmooth_markov_1, res_glmnet_undersmooth_markov_2, res_glmnet_undersmooth_markov_3)

truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds"))
sim_perf_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)

truth=truth[10,]

df <- res %>% filter(Target_parameter=="ATE", Estimator=="tmle") 
df$dataset_num <- 1:1500

df$cum_avg <- cummean(df$estimate)
df$cum_std.err <- cummean(df$std.err)


df <- df %>% mutate(bias=abs(cum_avg - truth$meanRD))
mean_bias <- abs(mean(df$estimate) - truth$meanRD)
head(df)

#df <- df %>% filter(dataset_num>1)

ggplot(df , aes(x=dataset_num, y=bias)) + geom_point() + geom_line() + 
  geom_hline(yintercept = mean_bias, linetype="dashed", color="red") + 
  geom_vline(xintercept = 1000, linetype="dashed") + 
  ggtitle("Mean Bias of ATE estimates as iterations increase") + xlab("Number of iterations") + ylab("Absolute bias")
#ggplot(df, aes(x=dataset_num, y=cum_std.err)) + geom_point() + geom_line() 

#Caption for this figure:
