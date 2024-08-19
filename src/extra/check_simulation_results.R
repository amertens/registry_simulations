shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)
library(ggplot2)
library(ggpubr)
u=lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
    do.call("require", list(X)) 
})
nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
truth<- readRDS(file=paste0(here::here(),"/data/sim_results/truth.rds"))
res = readRDS(file=paste0(here::here(),"/data/sim_results/sim_res.rds"))
setDT(res)
setDT(truth)
# restrict to analyses _2!? they have 500 repetitions. _1 have 50 and _3 have 200 ...
res2 <- res[grepl("_2",analysis)] %>% filter(!grepl("untruncated", estimator), grepl("markov", estimator))
# merge with truth
tt <- truth[time==10,.(Target_parameter=c("Risk(A=1)","Risk(A=0)","ATE","RelativeRisk"),
                       true_value=c(meanYa1,meanYa1-meanRD,meanRD,meanRR))]
res2 <- tt[res2,on="Target_parameter"]

# plot function to visualize the results
plot_simulation <- function(x,est){
    X <- x[Estimator=="tmle" & grepl(est,estimator)]
    xp=lapply(c("Risk(A=1)","Risk(A=0)","ATE","RelativeRisk"),
              function(tp){
                  ggplot(X[Target_parameter==tp],aes(y=estimate,x=estimator))+geom_boxplot()+geom_abline(intercept=X[Target_parameter==tp]$true_value[1],slope=0) + coord_cartesian(ylim=c(0,.02))
              })
    ggarrange(xp[[1]],xp[[2]],xp[[3]],xp[[4]],
              labels = c("Risk(A=1)","Risk(A=0)","ATE","RelativeRisk"),
              ncol = 2, nrow = 2)
}

plot_simulation(res2,est = "ridge")
# plot_simulation(res2,est = "glmnet")
 plot_simulation(res2,est = "EN")
