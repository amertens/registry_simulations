

calc_sim_performance <- function(res, truth, time=10){
  
  try(res=rbindlist(res))
  
  trueRR <- truth[time,2]
  trueRD <- truth[time,3]
  
  resRD <- res$estimate[res$Target_parameter=="ATE"]
  resRR <- res$estimate[res$Target_parameter=="RelativeRisk"]
  
  
  O_coverage_RD = mean(resRD-1.96*sd(resRD)< trueRD & trueRD < resRD+1.96*sd(resRD))  
  O_coverage_RR = mean(resRR-1.96*sd(resRR)< (trueRR) & trueRR < resRR+1.96*sd(resRR))  
  
  tab = data.frame(O_coverage_RD=O_coverage_RD,  O_coverage_RR=O_coverage_RR)
  
  return(tab)
}




calc_sim_performance <- function(WD=paste0(here::here(),"/sim_res/"), files, boot_iter_files=NULL, trueRR, trueRD, iptw=F){
  
  setwd(WD)
  d <- files %>% map(readRDS) %>% map_dfr(~bind_rows(.) , .id="analysis")
  d <- d %>% mutate(analysis = factor(analysis))
  levels(d$analysis) = files[as.numeric(levels(d$analysis))]
  d$analysis <- gsub(".RDS","",d$analysis)
  
  
  #transform iptw
  colnames(d)
  d.iptw <- d %>% select(analysis, starts_with("iptw.")) %>% mutate(analysis=paste0(analysis,"_iptw"))
  colnames(d.iptw) <- gsub("iptw.","",colnames(d.iptw))
  d.iptw <- d.iptw %>% rename(CI.2.5.=ci.lb, CI.97.5.=ci.ub) %>% mutate(IPTW="Yes")
  d <- d %>% select(!starts_with("iptw.")) %>% mutate(IPTW="No")
  
  summary(d$ate)
  summary(d.iptw$ate)
  
  if(iptw){
    d <- bind_rows(d, d.iptw)
  }
  
  d$true.RR <- trueRR
  d$true.RD <- trueRD
  
  # head(d)
  perf_tab_RR <- d %>% group_by(analysis) %>%
    mutate(variance=mean(((estimate)-mean((estimate)))^2),
           o.ci.lb = (estimate) - 1.96 * sqrt(variance),
           o.ci.ub = (estimate) + 1.96 * sqrt(variance)) %>%
    group_by(analysis) %>%
    summarize(
      bias=mean(log(estimate))-log(true.RR),
      variance=mean((mean(log(estimate))-log(estimate))^2),
      mse = bias^2 + variance,
      bias_se_ratio= bias/sqrt(variance),
      bias_se_ratio_emp= bias/mean(std.dev),
      coverage=mean(CI.2.5.<=true.RR & true.RR<=CI.97.5.)*100,
      #oracle coverage
      o.coverage=mean(o.ci.lb<=(true.RR) & (true.RR)<= o.ci.ub)*100,
      mean_ci_width=mean(log(CI.97.5.)-log(CI.2.5.)),
      power=mean((CI.2.5. > 1 & CI.97.5.>1)|(CI.2.5. < 1 & CI.97.5.<1))*100
    ) %>% filter(!is.na(variance)) %>%
    distinct()
  
  
  perf_tab_diff <- d %>% filter(!is.na(ate)) %>%
    subset(., select = -c(estimate, std.dev, CI.2.5., CI.97.5.)) %>%
    rename(estimate= ate, std.dev= ate.sd, CI.2.5.=ate.ci.lb, CI.97.5.=ate.ci.ub)%>%
    group_by(analysis) %>%
    mutate(variance=mean((estimate-mean(estimate))^2),
           o.ci.lb = estimate - 1.96 * sqrt(variance),
           o.ci.ub = estimate + 1.96 * sqrt(variance)) %>%
    mutate(#variance=mean((estimate-mean(estimate))^2),
      o.ci.lb = estimate - 1.96 * sd(estimate),
      o.ci.ub = estimate + 1.96 * sd(estimate)) %>%
    group_by(analysis) %>%
    summarize(
      mean_ate=mean(estimate),
      bias=mean((estimate))-(true.RD),
      variance=mean((estimate-mean(estimate))^2),
      mse = bias^2 + variance,
      bias_se_ratio= bias/sqrt(variance),
      bias_se_ratio_emp= bias/mean(std.dev),
      coverage=mean(CI.2.5.<=true.RD & true.RD<=CI.97.5.)*100,
      #oracle coverage
      o.coverage=mean(o.ci.lb<=true.RD & true.RD<= o.ci.ub)*100,
      mean_ci_width=mean((CI.97.5.)-(CI.2.5.)),
      power=mean((CI.2.5. > 0 & CI.97.5.>0)|(CI.2.5. < 0 & CI.97.5.<0))*100,
      o.power=mean((o.ci.lb > 0 & o.ci.ub>0)|(o.ci.lb < 0 & o.ci.ub<0))*100
    ) %>%
    distinct()
  
  if(!is.null(boot_iter_files)){
    setwd(paste0(here::here(),"/data/bootstrap/"))
    
    if( class(boot_iter_files) == "data.frame"){
      boot_iter_files$boot_iter <- as.character(1:nrow(boot_iter_files))
      boot_res <- boot_iter_files$boot_file %>% map(readRDS) %>% map_dfr(~bind_rows(.) , .id="boot_iter")
      boot_res <- left_join(boot_res, boot_iter_files, by=c("boot_iter"))
      boot_res <- boot_res %>% subset(., select = -c(boot_file))
    }else{
      boot_res <- boot_iter_files %>% map(readRDS) %>% map_dfr(~bind_rows(.) , .id="boot_iter")
      boot_res$analysis = "bootstrap"
    }
    
    
    
    #transform iptw
    if(iptw){
      boot_res.iptw <- boot_res %>% select(analysis, starts_with("iptw.")) %>% mutate(analysis=paste0(analysis,"_iptw"))
      colnames(boot_res.iptw) <- gsub("iptw.","",colnames(boot_res.iptw))
      boot_res.iptw <- boot_res.iptw %>% rename(CI.2.5.=ci.lb, CI.97.5.=ci.ub) %>% mutate(analysis=paste0(analysis,"-IPTW"))
      boot_res <- boot_res %>% select(!starts_with("iptw."))
      boot_res <- bind_rows(boot_res, boot_res.iptw)
    }
    
    
    if(!is.null(boot_res$analysis)){
      boot_CIs <- boot_res %>% group_by(iteration, analysis) %>%
        summarise(
          CI1=quantile(estimate,.025),
          CI2=quantile(estimate,.975),
          ate.CI1=quantile(ate,.025),
          ate.CI2=quantile(ate,.975)
        )
      
      perf_tab_RR_boot <- boot_CIs %>% group_by(analysis) %>%
        summarize(
          coverage = mean(CI1 <= d$true.RR[1] &  CI2 >= d$true.RR[1])*100,
          mean_ci_width=mean(log(CI2)-log(CI1)),
          power=mean((CI1 > 1 & CI2>1)|(CI1 < 1 & CI2<1))*100
        ) %>% as.data.frame()
      
      perf_tab_diff_boot <- boot_CIs %>% group_by(analysis) %>%
        summarize(
          coverage = mean(ate.CI1 <= d$true.RD[1] &  ate.CI2 >= d$true.RD[1])*100,
          mean_ci_width=mean((ate.CI2)-(ate.CI1)),
          power=mean((ate.CI1 > 0 & ate.CI2>0)|(ate.CI1 < 0 & ate.CI2<0))*100
        ) %>% as.data.frame()
      
    }else{
      boot_CIs <- boot_res %>% group_by(iteration) %>%
        summarise(
          CI1=quantile(estimate,.025),
          CI2=quantile(estimate,.975),
          ate.CI1=quantile(ate,.025),
          ate.CI2=quantile(ate,.975)
        )
      boot_CIs$analysis = "bootstrap"
      
      
      #calculate bootstrap performance
      perf_tab_RR_boot <- data.frame(
        coverage = mean(boot_CIs$CI1 <= d$true.RR[1] &  boot_CIs$CI2 >= d$true.RR[1])*100,
        mean_ci_width_logRR=mean(log(boot_CIs$CI2)-log(boot_CIs$CI1)),
        power=mean((boot_CIs$CI1 > 1 & boot_CIs$CI2>1)|(boot_CIs$CI1 < 1 & boot_CIs$CI2<1))*100
      )
      perf_tab_RR_boot
      
      perf_tab_diff_boot <- data.frame(
        coverage = mean(boot_CIs$ate.CI1 <= d$true.RD[1] &  boot_CIs$ate.CI2 >= d$true.RD[1])*100,
        mean_ci_width=mean((boot_CIs$ate.CI2)-(boot_CIs$ate.CI1)),
        power=mean((boot_CIs$ate.CI1 > 0 & boot_CIs$ate.CI2>0)|(boot_CIs$ate.CI1 < 0 & boot_CIs$ate.CI2<0))*100
      )
    }
    
    perf_tab_RR <- bind_rows(perf_tab_RR, perf_tab_RR_boot)
    perf_tab_diff <- bind_rows(perf_tab_diff, perf_tab_diff_boot)
  }
  
  #perf_tab_diff <- perf_tab_diff %>% filter(!grepl("iptw",analysis))
  
  perf_tab_RR <- clean_sim_results(perf_tab_RR)
  perf_tab_diff <- clean_sim_results(perf_tab_diff)
  
  perf_tab_RR$iptw <- ifelse(grepl("iptw", perf_tab_RR$filenames), "IPTW", "TMLE")
  perf_tab_diff$iptw <- ifelse(grepl("iptw", perf_tab_diff$filenames), "IPTW", "TMLE")
  
  return(list(perf_tab_RR=perf_tab_RR, perf_tab_diff=perf_tab_diff))
}



# # results cleaning
#
# d<-perf_tab_RR
# unique(d$analysis)

clean_sim_results <- function(d){
  
  #clean analysis
  d$filenames <- d$analysis
  d$analysis <- gsub("sim_res_","",d$analysis)
  d$analysis <- gsub("T4","",d$analysis)
  d$analysis <- gsub("T11","",d$analysis)
  
  
  #simulated data type
  d$simulated_data <- ifelse(grepl("outcome_blind_",d$analysis),"Outcome blind","Null")
  d$analysis <- gsub("outcome_blind_","",d$analysis)
  
  
  #variance estimator
  d <- d %>% mutate(variance_estimator = case_when(
    grepl("ic_",analysis) ~ "ic",
    grepl("tmle_",analysis) ~ "tmle",
    grepl("bootstrap_no_ties",analysis) ~ "bootstrap-no replacement",
    grepl("bootstrap_",analysis) ~ "bootstrap-clustered ID"))
  d$analysis <- gsub("ic_","",d$analysis)
  d$analysis <- gsub("tmle_","",d$analysis)
  d$analysis <- gsub("bootstrap","",d$analysis)
  
  
  #Qint
  d$Qint <- ifelse(grepl("Qint_",d$analysis),"Yes","No")
  d$analysis <- gsub("Qint_","",d$analysis)
  #DetQ
  d$DetQ <- ifelse(grepl("noDetQ_",d$analysis),"No","Yes")
  d$analysis <- gsub("noDetQ_","",d$analysis)
  #null_no_cens_sim_res_noDetQ_Qint_tmle_T4
  
  #censoring+competing risks
  d$censoring_in_data <- ifelse(grepl("no_cens_",d$analysis),"Uncensored","Censored+Competing risks")
  d$analysis <- gsub("cens_competing_risks_","",d$analysis)
  d$analysis <- gsub("no_cens_","",d$analysis)
  
  #Sl library
  d$estimator <- ifelse(grepl("glm",d$analysis),"GLM","LASSO")
  d$analysis <- gsub("glm_","",d$analysis)
  
  #sort dataframe
  d <- d %>% subset(., select = -c(analysis))
  #d <- d %>% select(analysis)
  
  d <- d[,which(colnames(d) %in% c("simulated_data", "estimator", "variance_estimator", "Qint","DetQ",
                                   "censoring_in_data","bias","variance","mse","bias_se_ratio","bias_se_ratio_emp",
                                   "coverage",        'o.coverage',  "abs.o.coverage",       'mean_ci_width',
                                   "power", "filenames", "IPTW"))]
  
  d <- d %>% arrange(o.coverage)
  
  
  return(d)
}


