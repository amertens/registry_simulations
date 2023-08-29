
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


#IC variance res
ic_res <- read.csv(file=paste0(here::here(),"/data/sim_perf_500reps.csv"))
ic_res=ic_res[ic_res$estimator=="ridge_undersmooth_markov",]
ic_res=ic_res %>% filter(Estimator=="tmle")

#TMLE variance res
truth<- readRDS(file=paste0(here::here(),"/data/sim_results/truth.rds"))
tmle_res <- readRDS(file=paste0(here::here(),"/data/sim_results/sim_res_undersmooth_ridge_markov_tmle.RDS"))

tmle_res = calc_sim_performance(
  res=tmle_res,
  truth=truth,
  time=10)
tmle_res=tmle_res %>% filter(Estimator=="tmle")

#Bootstrap variance res
boot1 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot1.RDS")) %>% mutate(boot_run=1)
boot2 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot2.RDS")) %>% mutate(boot_run=2)
boot3 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot3.RDS")) %>% mutate(boot_run=3)
boot4 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot4.RDS")) %>% mutate(boot_run=4)
boot5 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot5.RDS")) %>% mutate(boot_run=5)
boot_res <- bind_rows(boot1, boot2, boot3, boot4, boot5)
bootCIs <- boot_res %>% group_by(Target_parameter, sim_iter, boot_run) %>%
  summarise(
    boot.CI1=quantile(estimate,.025),
    boot.CI2=quantile(estimate,.975)
  )

boot_res_A1=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="Risk(A=1)"] < truth$meanYa1[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="Risk(A=1)"] > truth$meanYa1[10] )*100
boot_res_RD=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="ATE"] < truth$meanRD[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="ATE"] > truth$meanRD[10] )*100
boot_res_RR=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="RelativeRisk"] < truth$meanRR[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="RelativeRisk"] > truth$meanRR[10] )*100



res_tab <- data.frame(`Variance estimator`=c("Influence curve","Bootstrap","TMLE"),
                      `Y_{A=1} Coverage`=c(ic_res$coverage_Ya1, boot_res_A1,  tmle_res$coverage_Ya1),
                      `RD Coverage`=c(ic_res$coverage_RD, boot_res_RD, tmle_res$coverage_RD),
                      `RR Coverage`=c(ic_res$coverage_RR, boot_res_RR, tmle_res$coverage_RR))
knitr::kable(res_tab)

#Make two boxplots of bootstrap vs IC
bootCIs
