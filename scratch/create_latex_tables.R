

rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)
library(kableExtra)
library(tidyverse)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

truth<- readRDS(file=paste0(here::here(),"/data/sim_results/truth.rds"))
res = readRDS(file=paste0(here::here(),"/data/sim_results/sim_res.rds"))
res_null = readRDS(file=paste0(here::here(),"/data/sim_results/sim_res_null.rds"))

#calc performance
res_tab = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)
head(res_tab)

null_truth=truth
null_truth[10,3]=1
null_truth[10,4]=0

null_res_tab = calc_sim_performance(
  res=res,
  truth=null_truth,
  time=10)


#clean up labels
res_tab <- res_tab %>% subset(., select = -c(N_reps)) %>% rename(Algorithm=estimator,
                            `Y_{A=1} bias`=abs_bias_Ya1,
                            `Y_{A=1} variance`=mean_variance_Ya1,
                            `Y_{A=1} bias SE ratio`=bias_se_ratio_Ya1,
                            `RD bias`=abs_bias_RD,
                            `RD variance`=mean_variance_RD,
                            `RD bias SE ratio`=bias_se_ratio_RD,
                            `RR log-transformed bias`=abs_log_bias_RR,
                            `RR variance`=mean_variance_RR,
                            `RR bias SE ratio`=bias_se_ratio_RR,
                            `Y_{A=1} oracle 95% coverage`=O_coverage_Ya1,
                            `RD oracle 95% coverage`=O_coverage_RD,
                            `RR oracle 95% coverage`=O_coverage_RR)
null_res_tab <- null_res_tab %>% subset(., select = -c(N_reps)) %>% rename(Algorithm=estimator,
                                                               `Y_{A=1} bias`=abs_bias_Ya1,
                                                               `Y_{A=1} variance`=mean_variance_Ya1,
                                                               `Y_{A=1} bias SE ratio`=bias_se_ratio_Ya1,
                                                               `RD bias`=abs_bias_RD,
                                                               `RD variance`=mean_variance_RD,
                                                               `RD bias SE ratio`=bias_se_ratio_RD,
                                                               `RR log-transformed bias`=abs_log_bias_RR,
                                                               `RR variance`=mean_variance_RR,
                                                               `RR bias SE ratio`=bias_se_ratio_RR,
                                                               `Y_{A=1} oracle 95% coverage`=O_coverage_Ya1,
                                                               `RD oracle 95% coverage`=O_coverage_RD,
                                                               `RR oracle 95% coverage`=O_coverage_RR)





create_sim_latex_tab <- function(res_table, filename, measure="RR", bold=FALSE){
  
  res_table <- data.table(res_table)
  
  ## keep only markvov property estimates 
  res_table <- res_table[grep("markov",Algorithm),]
  
  ## adding rows for specifications of weights and lambda
  res_table[grep("untruncated",Algorithm),Truncation:="Untruncated"]
  res_table[is.na(Truncation),Truncation:="0.01"]
  
  res_table[grep("undersmooth",Algorithm),Lambda:="Undersmoothed"]
  res_table[is.na(Lambda) & grepl("ridge|glmnet|EN_",Algorithm)==TRUE,Lambda:="CV-minimum SE"]
  res_table[is.na(Lambda),Lambda:="N/A"]
  
  
  ##updating naming 
  res_table[,Estimator:=toupper(Estimator)]
  
  res_table[,alg_old := Algorithm]
  res_table[grep("glmnet",alg_old),Algorithm:="LASSO"]
  res_table[grep("EN_",alg_old),Algorithm:="Elastic Net"]
  res_table[grep("ridge",alg_old),Algorithm:="Ridge"]
  res_table[grep("glm_",alg_old),Algorithm:="GLM"]
    
  # sort 
  names(res_table)
  ##res_table <- setorder(res_table, -Estimator,-`RD oracle 95% coverage` ,`RD bias`)
  res_table <- setorder(res_table, -Estimator,-`RR oracle 95% coverage` ,`RR log-transformed bias`)

  head(res_table)
  # identify index of rows to highlight
  ## keep only estimates for the measure of interest (RR or RD)
  
                            if(measure=="RR"){
                              res_table <- res_table[,.(Estimator, 
                                                        Algorithm, 
                                                        Lambda, 
                                                        Truncation,                            
                                                        `RR oracle 95% coverage`,
                                                        `RR log-transformed bias`,
                                                        `RR variance`)]
                            }else if (measure=="RD"){
                              res_table <- res_table[,.(Estimator, 
                                                        Algorithm, 
                                                        Lambda, 
                                                        Truncation,
                                                        `RD oracle 95% coverage`,
                                                        `RD bias`,
                                                        `RD variance`)] 
                            }
  
  
  #save for html file
  res_xtable <- res_table %>%
    knitr::kable(
      format = "latex",
      align = "l",
      booktabs = TRUE,
      longtable = TRUE,
      linesep = "",
      digits =6) %>%
    kableExtra::kable_styling(
      position = "left"#,
      # latex_options = c("striped", "repeat_header"),
      # stripe_color = "gray!15"
    ) %>%
    kable_styling() 
  if(bold==TRUE){ 
    res_xtable <- res_xtable %>%
    row_spec(1, hline_after = T) %>%
    row_spec(1, bold=T,hline_after = T)
    } 
  
  
  save_kable(res_xtable, file=paste0(here::here(),"/tables/",filename,".tex"),float = FALSE,format="latex")
}

create_sim_latex_tab(null_res_tab, filename="RR_results_null_table", measure="RR",bold=FALSE)
create_sim_latex_tab(res_tab, filename="RR_results_sig_table", measure="RR",bold=TRUE)
create_sim_latex_tab(null_res_tab, filename="RD_results_null_table", measure="RD",bold=FALSE)
create_sim_latex_tab(res_tab, filename="RD_results_sig_table", measure="RD",bold=FALSE)

