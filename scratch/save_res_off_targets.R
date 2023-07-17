
glmnet_res_t1 <- tar_read(glmnet_res_t1)
glmnet_res_t2 <- tar_read(glmnet_res_t2)
glmnet_res_t4 <- tar_read(glmnet_res_t4)
glmnet_res_t6 <- tar_read(glmnet_res_t6)
glmnet_res_t8 <- tar_read(glmnet_res_t8)
glmnet_res <- tar_read(glmnet_res)
glmnet_undersmooth_res <- tar_read(glmnet_undersmooth_res)
glmnet_markov_res <- tar_read(glmnet_markov_res)
glmnet_markov_undersmooth_res <- tar_read(glmnet_markov_undersmooth_res)
ridge_res <- tar_read(ridge_res)
glm_res <- tar_read(glm_res)
glm_markov_res <- tar_read(glm_markov_res)

truth <- tar_read(truth)
truth$truth_df

glmnet_res_t1$analysis <- "glmnet_res_t1"
glmnet_res_t2$analysis <- "glmnet_res_t2"
glmnet_res_t4$analysis <- "glmnet_res_t4"
glmnet_res_t6$analysis <- "glmnet_res_t6"
glmnet_res_t8$analysis <- "glmnet_res_t8"
glmnet_res$analysis <- "glmnet_res"
glmnet_undersmooth_res$analysis <- "glmnet_undersmooth_res"
glmnet_markov_res$analysis <- "glmnet_markov_res"
glmnet_markov_undersmooth_res$analysis <- "glmnet_markov_undersmooth_res"
ridge_res$analysis <- "ridge_res"
glm_res$analysis <- "glm_res"
glm_markov_res$analysis <- "glm_markov_res"




save( glmnet_res_t1,
      glmnet_res_t2,
      glmnet_res_t4,
      glmnet_res_t6,
      glmnet_res_t8,
      glmnet_res,
      glmnet_undersmooth_res,
      glmnet_markov_res,
      glmnet_markov_undersmooth_res,
      ridge_res,
glm_res,
glm_markov_res, file=paste0(here::here(),"/scratch/sim_res_200_rep.Rdata"))



res <- bind_rows(
  glmnet_res,
  glmnet_undersmooth_res,
  glmnet_markov_res,
  glmnet_markov_undersmooth_res,
  ridge_res,
  glm_res,
  glm_markov_res) %>% filter(!is.na(estimate))


res_tmle <- res %>% filter(Target_parameter=="ATE", Estimator=="tmle")
res_tmle$true.RR=0.6940056  
res_tmle$true.RD= -0.005892922


perf_tab_diff <- res_tmle %>% 
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
    bias_se_ratio_emp= bias/mean(std.err),
    coverage=mean(lower<=true.RD & true.RD<=upper)*100,
    #oracle coverage
    o.coverage=mean(o.ci.lb<=true.RD & true.RD<= o.ci.ub)*100,
    mean_ci_width=mean((upper)-(lower)),
    power=mean((lower > 0 & upper>0)|(lower < 0 & upper<0))*100,
    o.power=mean((o.ci.lb > 0 & o.ci.ub>0)|(o.ci.lb < 0 & o.ci.ub<0))*100
  ) %>%
  distinct()
perf_tab_diff






res_tmle <- glmnet_res_t1 %>% filter(Target_parameter=="ATE", Estimator=="tmle")
res_tmle_RR <- glmnet_res_t1 %>% filter(Target_parameter=="RelativeRisk", Estimator=="tmle")

res_tmle$true.RR= 1.1226718  
res_tmle$true.RD= 0.000191


perf_tab_diff <- res_tmle %>% 
  group_by(analysis) %>%
  mutate(variance=mean((estimate-mean(estimate))^2),
         o.ci.lb = estimate - 1.96 * sqrt(variance),
         o.ci.ub = estimate + 1.96 * sqrt(variance)) %>%
  mutate(
    o.ci.lb = estimate - 1.96 * sd(estimate),
    o.ci.ub = estimate + 1.96 * sd(estimate)) %>%
  group_by(analysis) %>%
  summarize(
    mean_ate=mean(estimate),
    bias=mean((estimate))-(true.RD),
    variance=mean((estimate-mean(estimate))^2),
    mse = bias^2 + variance,
    bias_se_ratio= bias/sqrt(variance),
    bias_se_ratio_emp= bias/mean(std.err),
    coverage=mean(lower<=true.RD & true.RD<=upper)*100,
    o.coverage=mean(o.ci.lb<=true.RD & true.RD<= o.ci.ub)*100,
    mean_ci_width=mean((upper)-(lower)),
    power=mean((lower > 0 & upper>0)|(lower < 0 & upper<0))*100,
    o.power=mean((o.ci.lb > 0 & o.ci.ub>0)|(o.ci.lb < 0 & o.ci.ub<0))*100
  ) %>%
  distinct() %>% as.data.frame()
perf_tab_diff
