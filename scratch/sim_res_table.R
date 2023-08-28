

rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(knitr)
library(xtable)
library(kableExtra)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

nn=lapply(list.files("./reals/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)


ic_res <- read.csv(file=paste0(here::here(),"/data/sim_perf_500reps.csv"))
ic_res

#To do:
#add iptw
#remove non-markov

#temp
ic_res$Method <- "TMLE"

head(ic_res)
res_diff_coverage <- ic_res %>% filter(grepl("markov",estimator))
res_diff_coverage$filenames <- res_diff_coverage$estimator
res_diff_coverage$estimator <- gsub("_markov","",res_diff_coverage$estimator )
res_diff_coverage$truncation <- ifelse(grepl("untruncated",res_diff_coverage$estimator),"Untruncated","<0.01")
res_diff_coverage$penalty_selection <- ifelse(grepl("undersmooth",res_diff_coverage$estimator),"Undersmoothed","CV-minimum SE")
res_diff_coverage$estimator <- gsub("_untruncated","",res_diff_coverage$estimator )
res_diff_coverage$estimator <- gsub("EN","Elastic net",res_diff_coverage$estimator )
res_diff_coverage$estimator <- gsub("glmnet","Lasso",res_diff_coverage$estimator )
res_diff_coverage$estimator <- gsub("ridge","Ridge",res_diff_coverage$estimator )

res_table <- res_diff_coverage %>%
  subset(., select=c(Method, estimator, penalty_selection, truncation, abs_bias_RD, mean_variance_RD, bias_se_ratio_RD, O_coverage_RD)) %>%
  rename(Algorithm=estimator, 
         `g-bounds`=truncation,
         `Penalty selection`=penalty_selection,
         Bias=abs_bias_RD, 
         Variance=mean_variance_RD, 
         `Bias/SE`=bias_se_ratio_RD, 
         `Oracle coverage`=O_coverage_RD) %>%
  mutate(Bias=round(Bias,2),
         Variance=round(Variance,6),
         `Bias/SE`=round(`Bias/SE`,2)) %>%
  arrange(`Oracle coverage`)
 # %>% arrange(iptw,estimator,  Qint, DetQ,  bias, variance)


# identify index of rows to highlight
row.i.1 <- which(res_table$filenames=="glmnet_undersmooth_markov")

# 
# #save for html file
# res_table_protective <- res_table
# res_table_protective_raw <- res_diff_raw
# save(res_table_protective, res_table_protective_raw, file=paste0(here::here(),"/results/sim_performance_results_protective.Rdata"))

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


save_kable(res_xtable, file=paste0(here::here(),"/tables/sim_results_table.tex"),float = FALSE)

#Ask: should we add the oracle coverage for Ya1 and RR to the tables?

