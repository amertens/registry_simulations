
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


truth<- readRDS(file=paste0(here::here(),"/data/sim_results/truth.rds"))
ic_res <- read.csv(file=paste0(here::here(),"/data/sim_perf_500reps.csv"))


boot1 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot1.RDS")) %>% mutate(boot_run=1)
boot2 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot2.RDS")) %>% mutate(boot_run=2)
boot3 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot3.RDS")) %>% mutate(boot_run=3)
boot_res <- bind_rows(boot1, boot2, boot3)
bootCIs <- boot_res %>% group_by(Target_parameter, sim_iter, boot_run) %>%
  summarise(
    boot.CI1=quantile(estimate,.025),
    boot.CI2=quantile(estimate,.975)
  )

boot_res_A1=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="Risk(A=1)"] < truth$meanYa1[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="Risk(A=1)"] > truth$meanYa1[10] )*100
boot_res_RD=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="ATE"] < truth$meanRD[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="ATE"] > truth$meanRD[10] )*100
boot_res_RR=mean(bootCIs$boot.CI1[bootCIs$Target_parameter=="RelativeRisk"] < truth$meanRR[10] & bootCIs$boot.CI2[bootCIs$Target_parameter=="RelativeRisk"] > truth$meanRR[10] )*100

ic_res=ic_res[ic_res$estimator=="ridge_undersmooth_markov",]

res_tab <- data.frame(`Variance estimator`=c("Influence curve","Bootstrap","TMLE"),
                      `Y_{A=1} Coverage`=c(ic_res$coverage_Ya1, boot_res_A1, "100 (placeholder)"),
                      `RD Coverage`=c(ic_res$coverage_RD, boot_res_RD, 100),
                      `RR Coverage`=c(ic_res$coverage_RR, boot_res_RR, 100))
knitr::kable(res_tab)
