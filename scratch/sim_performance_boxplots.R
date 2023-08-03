


shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(data.table)
library(tidyverse)


# -------------------------------------------------------------------------------------------------------------
# load and clean results
# -------------------------------------------------------------------------------------------------------------


df <- bind_rows(
        data.frame(estimator="glm",tar_read( sim_res_tab_glm)),   
        data.frame(estimator="glmnet",tar_read( sim_res_tab_glmnet)),   
        data.frame(estimator="glmnet_undersmooth", tar_read(sim_res_tab_glmnet_undersmooth)),   
        #data.frame(tar_read(estimator="glmnet_markov", sim_res_tab_glmnet_markov)),   
        #data.frame(tar_read(estimator="glmnet_undersmooth_markov", sim_res_tab_glmnet_undersmooth_markov)),   
        data.frame(estimator="ridge", tar_read( sim_res_tab_ridge)),   
        data.frame(estimator="ridge_undersmooth", tar_read( sim_res_tab_ridge_undersmooth)),   
        #data.frame(tar_read(estimator="", sim_res_tab_ridge_markov)),   
        #data.frame(tar_read(estimator="", sim_res_tab_ridge_undersmooth_markov)),   
        data.frame(estimator="EN", tar_read(sim_res_tab_EN)),   
        data.frame(estimator="EN_undersmooth", tar_read(sim_res_tab_EN_undersmooth)),  
        #data.frame(tar_read(estimator="", sim_res_tab_EN_markov)), 
        #data.frame(tar_read(estimator="", sim_res_tab_EN_undersmooth_markov)),  

        data.frame(estimator="untruncated_g_glm", tar_read(sim_trunc_tab_glm)),   
        data.frame(estimator="untruncated_g_glmnet", tar_read(sim_trunc_tab_glmnet)),   
        data.frame(estimator="untruncated_g_glmnet_undersmooth",tar_read( sim_trunc_tab_glmnet_undersmooth)),   
        #data.frame(tar_read(estimator="", sim_trunc_tab_glmnet_markov)),   
        #data.frame(tar_read(estimator="", sim_trunc_tab_glmnet_undersmooth_markov)),  
        data.frame(estimator="untruncated_g_ridge",tar_read( sim_trunc_tab_ridge)),   
        data.frame(estimator="untruncated_g_ridge_undersmooth", tar_read(sim_trunc_tab_ridge_undersmooth)),   
        #data.frame(tar_read(estimator="", sim_trunc_tab_ridge_markov)),   
        #data.frame(tar_read(estimator="", sim_trunc_tab_ridge_undersmooth_markov)),   
        data.frame(estimator="untruncated_g_EN", tar_read(sim_trunc_tab_EN)),   
        data.frame(estimator="untruncated_g_EN_undersmooth", tar_read(sim_trunc_tab_EN_undersmooth))#,   
        #data.frame(tar_read(estimator="", sim_trunc_tab_EN_markov)),   
        #data.frame(tar_read(estimator="", sim_trunc_tab_EN_undersmooth_markov))
        )



df_null <- bind_rows(
                data.frame(tar_read(estimator="", null_res_tab_glm)),  
                data.frame(tar_read(estimator="", null_res_tab_glmnet)),   
                data.frame(tar_read(estimator="", null_res_tab_glmnet_undersmooth)),   
                data.frame(tar_read(estimator="", null_res_tab_glmnet_markov)),  
                data.frame(tar_read(estimator="", null_res_tab_glmnet_undersmooth_markov)),   
                data.frame(tar_read(estimator="", null_res_tab_ridge)),   
                data.frame(tar_read(estimator="", null_res_tab_ridge_undersmooth)),   
                data.frame(tar_read(estimator="", null_res_tab_ridge_markov)),   
                data.frame(tar_read(estimator="", null_res_tab_ridge_undersmooth_markov)),   
                data.frame(tar_read(estimator="", null_res_tab_EN)),   
                data.frame(tar_read(estimator="", null_res_tab_EN_undersmooth)),   
                data.frame(tar_read(estimator="", null_res_tab_EN_markov)), 
                data.frame(tar_read(estimator="", null_res_tab_EN_undersmooth_markov)))



# -------------------------------------------------------------------------------------------------------------
# Plot funcrions
# -------------------------------------------------------------------------------------------------------------



plot_bias <- function(res, xlimit=c(0.1, 2), target_parameter="ATE", truth=-0.0071862, Title="Estimator bias boxplots"){
  res <- res %>% filter(Target_parameter==target_parameter) %>% group_by(estimator) %>%
    mutate(mean_abs_bias=mean(abs(estimate-truth))) %>% ungroup() %>%
    arrange(mean_abs_bias) %>%
    mutate(estimator=factor(estimator, levels=unique(estimator)))
  
  
  p <- ggplot(res, aes(x=estimate, y=estimator, group=estimator)) + 
    geom_boxplot() + geom_vline(xintercept=truth, color="red") + 
    coord_cartesian(xlim=xlimit) + xlab("Estimate distribution (red line=truth") + ylab("Estimator (arranged by mean absolute bias)") +
    ggtitle(Title)
  
  if(target_parameter=="RelativeRisk"){
    p <- p + scale_x_continuous(trans = "log10")
  }
  
  return(p)
}

plot_se <- function(res, xlimit=c(0.1,0.5), target_parameter="ATE", Title="Estimator bias boxplots"){
  res <- res %>% filter(Target_parameter==target_parameter)
  res <- res %>% filter(Target_parameter==target_parameter) %>% group_by(estimator) %>%
    mutate(estimator_SE=sd(estimate)) %>% ungroup() %>%
    arrange(estimator_SE) %>%
    mutate(estimator=factor(estimator, levels=unique(estimator)))
  
  p <- ggplot(res, aes(x=std.err, y=estimator, group=estimator)) + 
    geom_boxplot() + 
    geom_point(aes(y=estimator, x=estimator_SE), color="red") + 
    coord_cartesian(xlim=xlimit) + 
    xlab("Observed SE distribution (red points=sd(estimate)") + ylab("Estimator (arranged by sd(estimate))") +
    ggtitle(Title)
  
  
  return(p)
}


# -------------------------------------------------------------------------------------------------------------
# Plot results
# -------------------------------------------------------------------------------------------------------------


bias_boxplots_ATE <- plot_bias(df, xlimit=c(-0.02, 0.005), Title="Estimator bias boxplots- Risk difference")
bias_boxplots_ATE
ggsave(bias_boxplots_ATE,  file=paste0(here::here(),"/figures/bias_boxplots_ATE.png"), height=5, width=8)

bias_boxplots_Ya1 <- plot_bias(df, target_parameter="Risk(A=1)", xlimit=c(0.001, 0.02), truth=0.0065,Title="Estimator bias boxplots- MeanY_(A=1)")
bias_boxplots_Ya1
ggsave(bias_boxplots_Ya1,  file=paste0(here::here(),"/figures/bias_boxplots_Ya1.png"), height=5, width=8)

bias_boxplots_RR <- plot_bias(df, target_parameter="RelativeRisk", xlimit=c(0.2,1.3), truth=0.5356286, Title="Estimator bias boxplots- Relative Risk")
bias_boxplots_RR
ggsave(bias_boxplots_RR,  file=paste0(here::here(),"/figures/bias_boxplots_RR.png"), height=5, width=8)

SE_boxplots_ATE <- plot_se(df, xlimit=c(0, 0.005), Title="Estimator bias boxplots- Risk difference")
SE_boxplots_ATE
ggsave(SE_boxplots_ATE,  file=paste0(here::here(),"/figures/SE_boxplots_ATE.png"), height=5, width=8)

SE_boxplots_Ya1 <- plot_se(df, xlimit=c(0, 0.005), target_parameter="Risk(A=1)", Title="Estimator bias boxplots- MeanY_(A=1)")
SE_boxplots_Ya1
ggsave(SE_boxplots_Ya1,  file=paste0(here::here(),"/figures/SE_boxplots_Ya1.png"), height=5, width=8)

SE_boxplots_RR <- plot_se(df, xlimit=c(0.1, 0.6), target_parameter="RelativeRisk", Title="Estimator bias boxplots- Relative Risk")
SE_boxplots_RR
ggsave(SE_boxplots_RR,  file=paste0(here::here(),"/figures/SE_boxplots_RR.png"), height=5, width=8)

