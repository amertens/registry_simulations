
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(data.table)
library(tidyverse)


#-------------------------------------------------------------
# SImulation performance
#-------------------------------------------------------------

sim_performance = tar_read(sim_performance)
knitr::kable(sim_performance)


ggplot(sim_performance, aes(x=mean_variance_RD , y=estimator_variance_RD )) + geom_point() + coord_cartesian(xlim=c(0,0.00001), ylim=c(0,0.00001)) +
  geom_abline(intercept = 0, slope = 1) +labs(x="Mean empirical variance",y="Estimator variance")
ggplot(sim_performance, aes(x=mean_variance_RR , y=estimator_variance_RR )) + geom_point() + geom_abline(intercept = 0, slope = 1) +labs(x="Mean empirical variance",y="Estimator variance")


#Make boxplot for each method, with 200 variances for each boxplot
#Show all 200


#-------------------------------------------------------------
# Plot of simulation performance convergence
#-------------------------------------------------------------

sim_res_glm_1000iter <- tar_read(sim_res_glm_1000iter)
sim_res_glm_1000iter_tab <- clean_sim_res(sim_res_glm_1000iter)
head(sim_res_glm_1000iter_tab)
tar_read(truth)

RR_res <- sim_res_glm_1000iter_tab %>% filter(Target_parameter=="RelativeRisk")
RR_res

res=tar_read(sim_res_tab_glmnet_undersmooth_markov)

plot_bias <- function(res, xlimit=c(0,1)){
  RR_res <- res %>% filter(Target_parameter=="RelativeRisk")
  RR_res
  p <- ggplot(RR_res, aes(x=estimate)) + geom_boxplot() + geom_vline(xintercept=0.5356286, color="red") + coord_cartesian(xlim=xlimit)
  return(p)
}

plot_se <- function(res, xlimit=c(0.1,0.5)){
  RR_res <- res %>% filter(Target_parameter=="RelativeRisk")
  RR_res
  p <- ggplot(RR_res, aes(x=std.err)) + geom_boxplot() + geom_vline(aes(xintercept=sd(estimate)), color="red") + coord_cartesian(xlim=xlimit)
  return(p)
}


sim_res_tab_glmnet=tar_read(sim_res_tab_glmnet)
plot_bias(sim_res_tab_glmnet)
plot_se(sim_res_tab_glmnet)

sim_res_tab_glmnet_undersmooth=tar_read(sim_res_tab_glmnet_undersmooth)
plot_bias(sim_res_tab_glmnet_undersmooth)
plot_se(sim_res_tab_glmnet_undersmooth)



RR_res <- res %>% filter(Target_parameter=="RelativeRisk")
RR_res
ggplot(RR_res, aes(x=estimate)) + geom_boxplot() + geom_vline(xintercept=0.5356286, color="red")
ggplot(RR_res, aes(x=std.err)) + geom_boxplot() + geom_vline(aes(xintercept=sd(estimate)), color="red")

sd(RR_res$estimate)

#in the register data, using glmnet, no reason to fear that the CI are anti-conservative even though we have rare outcome/


ATE_res <- sim_res_glm_1000iter_tab %>% filter(Target_parameter=="ATE")


calc_mean_bias <- function(res, nsub=100, truth= -0.0071862){
  res_sub = res[sample(1:nrow(res),replace=FALSE,size=nsub)]
            bias=mean(abs(res_sub$estimate-truth))
            mean_variance_RD=mean((res_sub$std.err)^2)
            estimator_variance_RD=mean(((res_sub$estimate)-mean((res_sub$estimate)))^2)
            coverage_RD=mean(res_sub$lower<=truth & truth<=res_sub$upper)*100
            O_coverage_RD=mean(res_sub$estimate-1.96*sd(res_sub$estimate)< truth & truth < res_sub$estimate+1.96*sd(res_sub$estimate))*100
  

  return(data.frame(bias=bias,mean_variance_RD=mean_variance_RD,estimator_variance_RD=estimator_variance_RD,coverage_RD=coverage_RD,O_coverage_RD=O_coverage_RD))
}

calc_mean_bias(res=ATE_res, nsub=100, truth= -0.0071862)

bias_df=rbindlist(lapply(1:1000, function(x) calc_mean_bias(res=ATE_res, nsub=x, truth= -0.0071862)))
bias_df$n <- 1:1000

ggplot(bias_df, aes(x=n, y=bias)) + geom_line() + geom_smooth(se=F)
ggplot(bias_df, aes(x=n, y=mean_variance_RD)) + geom_line() + geom_smooth(se=F)
ggplot(bias_df, aes(x=n, y=estimator_variance_RD)) + geom_line() + geom_smooth(se=F)
ggplot(bias_df, aes(x=n, y=O_coverage_RD)) + geom_line() + geom_smooth(se=F) + coord_cartesian(ylim=c(90,100)) + geom_hline(yintercept = 95, linetype="dashed",color="red")
ggplot(bias_df, aes(x=n, y=coverage_RD )) + geom_line() + geom_smooth(se=F) + coord_cartesian(ylim=c(80,100)) + geom_hline(yintercept = 95, linetype="dashed",color="red")

#Make this with all the actual results
#ggplot(bias_df, aes(x=mean_variance_RD, y=estimator_variance_RD)) + geom_point(alpha=0.5) + geom_smooth(se=F) + geom_abline(intercept = 0, slope = 1)

#-------------------------------------------------------------
# Compare IC and bootstrap CI width
#-------------------------------------------------------------

sim_res_glm_bootstrap = tar_read(sim_res_glm_bootstrap)

boot_CIs <- sim_res_glm_bootstrap %>% group_by(tar_batch , Target_parameter) %>%
  summarise(
    CI1=quantile(estimate,.025),
    CI2=quantile(estimate,.975)
  ) %>% filter(Target_parameter == "ATE")

mean(boot_CIs$CI2 -boot_CIs$CI1)


sim_res_tab_glm = tar_read(sim_res_tab_glm) %>% filter(Target_parameter == "ATE")
mean(sim_res_tab_glm$upper -sim_res_tab_glm$lower)

mean(boot_CIs$CI2 > (-0.0071862) & (-0.0071862) > boot_CIs$CI1)
mean(sim_res_tab_glm$upper > (-0.0071862) & (-0.0071862) > sim_res_tab_glm$lower)



#-------------------------------------------------------------
# Compare ridge and lasso estimates
#-------------------------------------------------------------

sim_res_tab_glmnet_undersmooth_markov = tar_read(sim_res_tab_glmnet_undersmooth_markov)  %>% filter(Target_parameter == "ATE") %>% arrange(estimate)
sim_res_tab_ridge_undersmooth_markov = tar_read(sim_res_tab_ridge_undersmooth_markov)  %>% filter(Target_parameter == "ATE") %>% arrange(estimate)

plotdf = data.frame(x=sim_res_tab_glmnet_undersmooth_markov$estimate, y=sim_res_tab_ridge_undersmooth_markov$estimate)
ggplot(plotdf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept = 0, slope = 1) 


#-------------------------------------------------------------
# Plot of simulation performance convergence
#-------------------------------------------------------------

tar_comp_lasso1 = clean_sim_res(tar_read(tar_comp_lasso1))
tar_comp_ridge1 = clean_sim_res(tar_read(tar_comp_ridge1))
tar_comp_lasso2 = clean_sim_res(tar_read(tar_comp_lasso2))
tar_comp_ridge2 = clean_sim_res(tar_read(tar_comp_ridge2))
tar_comp_lasso3 = clean_sim_res(tar_read(tar_comp_lasso3))
tar_comp_ridge3 = clean_sim_res(tar_read(tar_comp_ridge3))
tar_comp_lasso4 = clean_sim_res(tar_read(tar_comp_lasso4))
tar_comp_ridge4 = clean_sim_res(tar_read(tar_comp_ridge4))
tar_comp_lasso5 = clean_sim_res(tar_read(tar_comp_lasso5))
tar_comp_ridge5 = clean_sim_res(tar_read(tar_comp_ridge5))


df <- bind_rows(
  merge(tar_comp_lasso1, tar_comp_ridge1, by = c("Target_parameter", "Estimator", "iteration")),
  merge(tar_comp_lasso2, tar_comp_ridge2, by = c("Target_parameter", "Estimator", "iteration")),
  merge(tar_comp_lasso3, tar_comp_ridge3, by = c("Target_parameter", "Estimator", "iteration")),
  merge(tar_comp_lasso4, tar_comp_ridge4, by = c("Target_parameter", "Estimator", "iteration")),
  merge(tar_comp_lasso5, tar_comp_ridge5, by = c("Target_parameter", "Estimator", "iteration"))
  )

df_ate <- df %>% filter(Target_parameter =="ATE") %>% mutate(diff=estimate.x - estimate.y, pval_change=(pvalue.x-pvalue.y)/pvalue.x*100) 
df_RR <- df %>% filter(Target_parameter =="RelativeRisk") %>% mutate(diff=estimate.x - estimate.y, pval_change=(pvalue.x-pvalue.y)/pvalue.x*100)

ggplot(df_ate, aes(x=estimate.x, y=estimate.y)) + geom_point() + geom_abline(intercept = 0, slope = 1) +labs(x="lasso",y="ridge")
ggplot(df_RR, aes(x=estimate.x, y=estimate.y)) + geom_point()+labs(x="lasso",y="ridge")
