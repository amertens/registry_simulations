

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



#Need to updated and apply the below code (doesn't need to be a function)
create_sim_latex_tab <- function(res_table){
  
  # identify index of rows to highlight
  row.i.1 <- which(res_table$filenames=="sim_res_DetQ_ic_v3")
  
  res_table <- res_table %>% select(filenames,  iptw,
                                    estimator, Qint,
                                    bias,variance,bias_se_ratio,oracle.coverage) %>%
    rename(Algorithm=estimator,  `Q-int`=Qint,  Estimator=iptw, `Bias/SE`=bias_se_ratio,
           Bias=bias, Variance=variance, `Oracle coverage`=oracle.coverage)
  
  row.i.1 <- which(res_table$filenames=="sim_res_DetQ_ic_v3")
  
  
  #save for html file
  res_table_protective <- res_table
  res_table_protective_raw <- res_diff_raw
  
  res_table <- res_table %>% subset(., select = -c(filenames))
  
  print(as.data.frame(res_table))
  
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
    kable_styling()%>%
    row_spec(row.i.1-1, hline_after = T) %>%
    row_spec(row.i.1, bold=T,hline_after = T)
  
  
  #save_kable(res_xtable, file="C:/Users/andre/Documents/jici/diab-dementia-server-code/tables/sim_results_table.tex",float = FALSE)
  
}
