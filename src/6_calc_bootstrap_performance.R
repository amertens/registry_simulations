
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)
library(kableExtra)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)


#IC variance res
ic_res <- read.csv(file=paste0(here::here(),"/data/sim_perf_1000reps.csv"))
ic_res=ic_res[ic_res$estimator=="res_ridge_undersmooth",]
ic_res=ic_res %>% filter(Estimator=="tmle")

#TMLE variance res
truth=readRDS(paste0(here::here(),"/data/sim_results/truth.rds")) %>% subset(., select=-c(meanYa0))

tmle_res <- readRDS(file=paste0(here::here(),"/data/sim_results/sim_res_tmle_var.RDS"))

tmle_res %>% group_by(Target_parameter,Estimator) %>% summarize(mean(estimate))


tmle_res = calc_sim_performance(
  res=tmle_res,
  truth=truth,
  time=10)
tmle_res=tmle_res %>% filter(Estimator=="tmle")


#Bootstrap variance res first 500 reps
boot1 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot1.RDS")) %>% mutate(boot_run=1)
boot2 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot2.RDS")) %>% mutate(boot_run=2)
boot3 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot3.RDS")) %>% mutate(boot_run=3)
boot4 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot4.RDS")) %>% mutate(boot_run=4)
boot5 <- readRDS(file=paste0(here::here(),"/data/sim_results/res_ridge_undersmooth_markov_boot5.RDS")) %>% mutate(boot_run=5)
boot_res <- bind_rows(boot1, boot2, boot3, boot4, boot5)
bootCIs <- boot_res %>% group_by(Target_parameter, sim_iter, boot_run) %>%
  summarise(
    std.err=sd(estimate),
    boot.CI1=quantile(estimate,.025),
    boot.CI2=quantile(estimate,.975)
  )

bootCI_ATE = bootCIs %>% filter(Target_parameter=="ATE")  %>% arrange(boot.CI2+boot.CI1)
bootCI_ATE$rep =1:500
ggplot(bootCI_ATE, aes(x=rep)) +
  geom_hline(yintercept=(-0.0071862)) +
  geom_linerange(aes(ymin=boot.CI1, ymax=boot.CI2),alpha=0.5) 

mean(bootCI_ATE$boot.CI1< truth$meanRD[10] & bootCI_ATE$boot.CI2> truth$meanRD[10] )*100
summary(bootCI_ATE$boot.CI1)
summary(bootCI_ATE$boot.CI2)

#boot results -new 500 reps
bootCIs2 <- readRDS(file=paste0(here::here(),"/data/sim_results/bootCIs_501_1000.rds")) %>% mutate(boot_run=9, dataset_num=dataset_num+500)


bootCI_df <- bind_rows(bootCIs, bootCIs2) 
dim(bootCI_df %>% filter(Target_parameter=="ATE"))

boot_res_A1=mean(bootCI_df$boot.CI1[bootCI_df$Target_parameter=="Risk(A=1)"] < truth$meanYa1[10] & bootCI_df$boot.CI2[bootCI_df$Target_parameter=="Risk(A=1)"] > truth$meanYa1[10] )*100
boot_res_RD=mean(bootCI_df$boot.CI1[bootCI_df$Target_parameter=="ATE"] < truth$meanRD[10] & bootCI_df$boot.CI2[bootCI_df$Target_parameter=="ATE"] > truth$meanRD[10] )*100
boot_res_RR=mean(bootCI_df$boot.CI1[bootCI_df$Target_parameter=="RelativeRisk"] < truth$meanRR[10] & bootCI_df$boot.CI2[bootCI_df$Target_parameter=="RelativeRisk"] > truth$meanRR[10] )*100

res_tab <- data.frame(`Variance estimator`=c("Influence curve","Bootstrap","TMLE"),
                      `Y_{A=1} Coverage`=c(ic_res$coverage_Ya1, boot_res_A1, tmle_res$coverage_Ya1),
                      `RD Coverage`=c(ic_res$coverage_RD, boot_res_RD, tmle_res$coverage_RD),
                      `RR Coverage`=c(ic_res$coverage_RR, boot_res_RR, tmle_res$coverage_RR))
colnames(res_tab) <- c("Variance estimator","Y_{A=1} Coverage","RD Coverage","RR Coverage")
knitr::kable(res_tab, format="simple")



res_xtable <- res_tab %>%
  dplyr::select(`Variance estimator`, `RD Coverage`) %>%
  knitr::kable(
    format = "latex",
    align = "l",
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    digits =6) %>%
  kableExtra::kable_styling(
    position = "left"#
  )

save_kable(res_xtable, file=paste0(here::here(),"/tables/variance_comp_tab.tex"),float = FALSE,format="latex")



#show undersmoothed ridge versus CV lasso
