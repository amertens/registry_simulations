shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)
u=lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})
nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
truth<- readRDS(file=paste0(here::here(),"/data/sim_results/truth.rds"))
res = readRDS(file=paste0(here::here(),"/data/sim_results/sim_res.rds"))
setDT(res)
setDT(truth)
res2 <- res[grepl("_2",analysis)]
ridge <- res2[Estimator=="tmle" & grepl("ridge",estimator)]
tt <- truth[time==10,.(Target_parameter=c("Risk(A=1)","Risk(A=0)","ATE","RelativeRisk"),
                       true_value=c(meanYa1,meanYa1-meanRD,meanRD,meanRR))]
ridge <- tt[ridge,on="Target_parameter"]
## ridge[,table(analysis)]
xp=lapply(tt$Target_parameter,function(tp){
  ggplot(ridge[Target_parameter==tp],aes(y=estimate,x=estimator))+geom_boxplot()+geom_abline(intercept=ridge[Target_parameter==tp]$true_value[1],slope=0)
})
# Boxplot for Risk(A=1) vs truth
xp[[1]]
# Boxplot for Risk(A=0) vs truth
xp[[2]]
# Boxplot for ATE vs truth
xp[[3]]
# Boxplot for RR vs truth
xp[[4]]
