Diabetes dementia analysis simulation results
================
Andrew Mertens

## Simulation setup

- N: 100,000

- Iterations: 500

- True $Y_{\bar{A}=1}$: 0.0082798

- True RD: -0.0071862

- True RR: 0.5356286

## Calculate simulation performance

``` r
res = readRDS(file=paste0(here::here(),"/data/sim_results/sim_res.rds"))

ic_res = calc_sim_performance(
  res=res,
  truth=truth,
  time=10)
```

    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.
    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.
    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.

``` r
head(ic_res)
```

    ##   Estimator                         estimator N_reps abs_bias_Ya1
    ## 1      iptw                                EN    500  0.003229815
    ## 2      iptw                         EN_markov    500  0.003229815
    ## 3      iptw             EN_markov_untruncated    500  0.003229815
    ## 4      iptw                    EN_undersmooth    500  0.003229815
    ## 5      iptw             EN_undersmooth_markov    500  0.003229815
    ## 6      iptw EN_undersmooth_markov_untruncated    500  0.003123814
    ##   estimator_variance_Ya1 mean_variance_Ya1 bias_se_ratio_Ya1 coverage_Ya1
    ## 1           5.025682e-06      2.429615e-06          2.072091         48.4
    ## 2           5.025682e-06      2.308602e-06          2.125705         47.4
    ## 3           5.025682e-06      2.413355e-06          2.079060         48.2
    ## 4           5.025682e-06      4.127306e-06          1.589807         67.2
    ## 5           5.025682e-06      4.052985e-06          1.604317         66.6
    ## 6           5.098092e-06      9.014282e-06          1.040446         80.6
    ##   O_coverage_Ya1 abs_bias_RD estimator_variance_RD mean_variance_RD
    ## 1           75.0 0.004278272          5.390669e-06     2.738929e-06
    ## 2           75.0 0.004278272          5.390669e-06     2.615863e-06
    ## 3           75.0 0.004278272          5.390669e-06     2.720630e-06
    ## 4           75.0 0.004278272          5.390669e-06     4.246168e-06
    ## 5           75.0 0.004278272          5.390669e-06     4.147353e-06
    ## 6           75.8 0.004422596          5.564185e-06     9.108648e-06
    ##   bias_se_ratio_RD coverage_RD O_coverage_RD abs_log_bias_RR
    ## 1         2.585104        35.0            54       0.1732900
    ## 2         2.645214        34.6            54       0.1732900
    ## 3         2.593783        35.2            54       0.1732900
    ## 4         2.076203        44.6            54       0.1732900
    ## 5         2.100791        43.4            54       0.1732900
    ## 6         1.465380        57.0            52       0.1822519
    ##   estimator_variance_RR mean_variance_RR bias_se_ratio_RR coverage_RR
    ## 1           0.009733954       0.02001366        1.2249270        79.6
    ## 2           0.009733954       0.01904607        1.2556562        78.2
    ## 3           0.009733954       0.01984408        1.2301497        79.2
    ## 4           0.009733954       0.03332591        0.9492540        89.8
    ## 5           0.009733954       0.03260497        0.9596913        89.2
    ## 6           0.010011311       0.07178766        0.6802166        92.2
    ##   O_coverage_RR
    ## 1          93.8
    ## 2          93.8
    ## 3          93.8
    ## 4          93.8
    ## 5          93.8
    ## 6          93.6

## Simulation results

#### Comparison of IPTW and TMLE estimators

Both from Lasso models with $\lambda$ selected at the minimum
cross-validated SE. All tables below are sorted by the RD bias.

| Estimator | Algorithm | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | glmnet    |       0.00188 |                 0 |                1.23266 |                         85.6 | 0.00254 |           0 |          1.56835 |                   73.2 |                 0.31340 |     0.05211 |          1.37291 |                   82.6 |
| iptw      | glmnet    |       0.00323 |                 0 |                2.09748 |                         75.0 | 0.00428 |           0 |          2.61471 |                   54.0 |                 0.17329 |     0.01954 |          1.23971 |                   93.8 |

Result: TMLE performs better

#### Comparison of algorithms

All penalized regressions with $\lambda$ selected with undersmoothing,
truncation of $g$ at \< 0.01, and using the t-1 prior L nodes in $g$ and
$Q$ estimation.

RF= random forest

| Estimator | Algorithm                 | Algorithm Alpha | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:--------------------------|:----------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | ridge_undersmooth_markov  | Ridge           |       0.00182 |                 0 |                0.82860 |                         93.2 | 0.00185 |           0 |          0.83441 |                   95.2 |                 0.20153 |     0.06036 |          0.82026 |                   96.2 |
| tmle      | EN_undersmooth_markov     | Elastic Net     |       0.00196 |                 0 |                0.99201 |                         93.0 | 0.00196 |           0 |          0.97946 |                   95.6 |                 0.21502 |     0.04786 |          0.98285 |                   95.2 |
| tmle      | glmnet_undersmooth_markov | Lasso           |       0.00196 |                 0 |                0.99170 |                         93.0 | 0.00196 |           0 |          0.97899 |                   95.6 |                 0.21512 |     0.04796 |          0.98233 |                   95.2 |
| tmle      | glm                       | NA              |       0.00214 |                 0 |                1.05802 |                         92.0 | 0.00204 |           0 |          0.97350 |                   95.0 |                 0.22221 |     0.04837 |          1.01038 |                   95.8 |
| tmle      | RF                        | NA              |       0.00319 |                 0 |                2.96372 |                         86.4 | 0.00277 |           0 |          2.32506 |                   92.0 |                 0.28096 |     0.01065 |          2.72232 |                   92.2 |

Result: Ridge performs best

#### Comparison of penalized regression setup

All from Ridge models

<!-- = $\lambda$ selected at the minimum cross-validated SE -->
<!-- penalized regressions with $\lambda$ selected at the minimum cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L nodes in $g$ and $Q$ estimation. -->

| Undersmoothed        | Truncation  | Markov      | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:---------------------|:------------|:------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| Undersmoothed Lambda | g \< 0.01   | Markov L    |       0.00182 |             0e+00 |                0.82860 |                         93.2 | 0.00185 |       0e+00 |          0.83441 |                   95.2 |                 0.20153 |     0.06036 |          0.82026 |                   96.2 |
| Undersmoothed Lambda | g \< 0.01   | All prior L |       0.00177 |             0e+00 |                0.80592 |                         93.8 | 0.00187 |       0e+00 |          0.84377 |                   94.4 |                 0.20297 |     0.06262 |          0.81115 |                   95.6 |
| Undersmoothed Lambda | Untruncated | Markov L    |       0.00224 |             1e-05 |                0.68394 |                         93.6 | 0.00224 |       1e-05 |          0.68329 |                   96.2 |                 0.24283 |     0.10911 |          0.73515 |                   94.8 |
| Undersmoothed Lambda | Untruncated | All prior L |       0.00228 |             1e-05 |                0.69416 |                         93.8 | 0.00233 |       1e-05 |          0.70866 |                   95.0 |                 0.25208 |     0.11001 |          0.76003 |                   95.0 |
| Min SE Lambda        | Untruncated | Markov L    |       0.00196 |             0e+00 |                0.96693 |                         87.4 | 0.00255 |       0e+00 |          1.20798 |                   76.0 |                 0.32162 |     0.08782 |          1.08529 |                   83.2 |
| Min SE Lambda        | Untruncated | All prior L |       0.00197 |             0e+00 |                0.94833 |                         87.8 | 0.00258 |       0e+00 |          1.19462 |                   77.2 |                 0.32306 |     0.09047 |          1.07408 |                   83.2 |
| Min SE Lambda        | g \< 0.01   | Markov L    |       0.00197 |             0e+00 |                1.07467 |                         81.6 | 0.00261 |       0e+00 |          1.35587 |                   68.6 |                 0.33011 |     0.07937 |          1.17176 |                   79.4 |
| Min SE Lambda        | g \< 0.01   | All prior L |       0.00197 |             0e+00 |                1.07160 |                         81.8 | 0.00265 |       0e+00 |          1.36925 |                   66.4 |                 0.33284 |     0.08024 |          1.17497 |                   79.4 |

Result: Undersmoothed Ridge with default truncation and markov process
for L works best

## Comparison of IC and Bootstrap based coverages

    ## `summarise()` has grouped output by 'Target_parameter', 'sim_iter'. You can
    ## override using the `.groups` argument.

| Variance.estimator | Y\_.A.1..Coverage | RD.Coverage | RR.Coverage |
|:-------------------|------------------:|------------:|------------:|
| Influence curve    |              94.2 |        89.2 |        94.2 |
| Bootstrap          |              94.2 |        94.8 |        95.4 |

Results: Bootstrap has better coverage for RR and especially RD, but
slightly worse that for the mean of $Y_{A=1}$

<!-- ## Plots of IC vs Bootstrap coverage -->
<!-- For the coverage of the risk difference (Red line is true RD) -->
<!-- ```{r} -->
<!--  ic_plot_df <- res %>% filter(Estimator=="tmle",Target_parameter=="ATE",grepl("ridge_undersmooth_markov", estimator),!grepl("untrunc", estimator)) -->
<!--  ggplot(ic_plot_df %>% arrange(lower) %>% -->
<!--           mutate(iteration=row_number()), aes(x=iteration)) + -->
<!--           geom_linerange(aes(ymin=lower, ymax=upper), alpha=0.5) + -->
<!--           coord_flip() + -->
<!--           geom_hline(yintercept = truth[10,4]) + ggtitle("IC") -->
<!--  ggplot(bootCIs %>% filter(Target_parameter=="ATE") %>% arrange(boot.CI1) %>% ungroup() %>% -->
<!--           mutate(iteration=row_number()), aes(x=iteration)) + -->
<!--          geom_linerange(aes(ymin=boot.CI1, ymax=boot.CI2), alpha=0.5) + -->
<!--          coord_flip() + -->
<!--          geom_hline(yintercept = truth[10,4]) + ggtitle("Bootstrap") -->
<!-- ``` -->
<!-- (Notes, feel free to give visualization tips for comparing the bootstrap and IC. Like should I do boxplots of the SE's?) -->

## Updated needed in the registry analysis

Because of updated in (I think) TMLE_for_breakfast, there are now D
nodes for competing events, and when we run the old code we get the
error:

Current error: “Cannot both specify determistic.Q.nodes and Dnodes”

I think I just need to update both the run_tmle() and prepare_tmle()
functions to be able to pass a Dnodes=NULL argument through the
functions, because right now Dnodes are automatically set to be any
variable with “dead” within the name. But do I do that on Statistics
Denmark or make the changes on Ltmle_for_breakfast, and if the latter,
how do I export the update so it’s on Statistics Denmark?

## Estimator to run:

undersmoothed ridge with markov=TRUE

``` r
v=run_ltmle(OUT="dementia",k=10,test=FALSE,treatment_regimens=treatment_regimens,outcomes=outcomes,baseline_covariates=baseline_covariates,timevar_covariates=timevar_covariates,det.Q.function=det.Q.function,abar="try",
            SL.library="glmnet",SL.cvControl=list(selector="undersmooth",alpha=1), 
            Markov=Markov_vars)
```

With Markov_vars being all time-varying L nodes
