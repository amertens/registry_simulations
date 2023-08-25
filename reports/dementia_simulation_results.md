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
    ## 1      iptw                                EN    250  0.003147734
    ## 2      iptw                         EN_markov    250  0.003147734
    ## 3      iptw             EN_markov_untruncated    250  0.003147734
    ## 4      iptw                    EN_undersmooth    250  0.003147734
    ## 5      iptw             EN_undersmooth_markov    250  0.003147734
    ## 6      iptw EN_undersmooth_markov_untruncated    226  0.003124008
    ##   estimator_variance_Ya1 mean_variance_Ya1 bias_se_ratio_Ya1 coverage_Ya1
    ## 1           4.716835e-06      2.409565e-06          2.027816     48.00000
    ## 2           4.716835e-06      2.292094e-06          2.079130     46.40000
    ## 3           4.716835e-06      2.373743e-06          2.043060     47.20000
    ## 4           4.716835e-06      4.098073e-06          1.554920     66.40000
    ## 5           4.716835e-06      4.015761e-06          1.570775     66.80000
    ## 6           4.825249e-06      9.723021e-06          1.001870     83.18584
    ##   O_coverage_Ya1 abs_bias_RD estimator_variance_RD mean_variance_RD
    ## 1       73.60000 0.004374549          5.062876e-06     2.719256e-06
    ## 2       73.60000 0.004374549          5.062876e-06     2.599669e-06
    ## 3       73.60000 0.004374549          5.062876e-06     2.681334e-06
    ## 4       73.60000 0.004374549          5.062876e-06     4.225610e-06
    ## 5       73.60000 0.004374549          5.062876e-06     4.110316e-06
    ## 6       76.54867 0.004422784          5.098313e-06     9.817649e-06
    ##   bias_se_ratio_RD coverage_RD O_coverage_RD abs_log_bias_RR
    ## 1         2.652823    31.60000      51.60000       0.1692712
    ## 2         2.713153    30.80000      51.60000       0.1692712
    ## 3         2.671516    31.60000      51.60000       0.1692712
    ## 4         2.128083    42.80000      51.60000       0.1692712
    ## 5         2.157723    41.60000      51.60000       0.1692712
    ## 6         1.411536    57.52212      51.32743       0.1723349
    ##   estimator_variance_RR mean_variance_RR bias_se_ratio_RR coverage_RR
    ## 1           0.009055339       0.02000336        1.1968277    77.60000
    ## 2           0.009055339       0.01906596        1.2258965    76.40000
    ## 3           0.009055339       0.01965929        1.2072557    77.20000
    ## 4           0.009055339       0.03356044        0.9239943    92.00000
    ## 5           0.009055339       0.03278173        0.9349045    90.40000
    ## 6           0.009189393       0.07611912        0.6246349    93.36283
    ##   O_coverage_RR
    ## 1      92.40000
    ## 2      92.40000
    ## 3      92.40000
    ## 4      92.40000
    ## 5      92.40000
    ## 6      92.47788

## Simulation results

#### Comparison of IPTW and TMLE estimators

Both from Lasso models with $\lambda$ selected at the minimum
cross-validated SE. All tables below are sorted by the RD bias.

| Estimator | Algorithm | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | glmnet    |       0.00189 |                 0 |                1.23854 |                     86.26667 | 0.00255 |           0 |          1.57448 |               71.33333 |                 0.31399 |     0.05204 |          1.37645 |                   82.0 |
| iptw      | glmnet    |       0.00315 |                 0 |                2.05372 |                     73.60000 | 0.00437 |           0 |          2.68440 |               51.60000 |                 0.16927 |     0.01950 |          1.21220 |                   92.4 |

Result: TMLE performs better

#### Comparison of algorithms

All penalized regressions with $\lambda$ selected with undersmoothing,
truncation of $g$ at \< 0.01, and using the t-1 prior L nodes in $g$ and
$Q$ estimation.

RF= random forest

| Estimator | Algorithm                 | Algorithm Alpha | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:--------------------------|:----------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | ridge_undersmooth_markov  | Ridge           |       0.00182 |                 0 |                0.82617 |                     93.46667 | 0.00184 |           0 |          0.82992 |               95.20000 |                 0.20108 |     0.06054 |          0.81722 |               95.73333 |
| tmle      | EN_undersmooth_markov     | Elastic Net     |       0.00197 |                 0 |                0.99264 |                     92.93333 | 0.00195 |           0 |          0.97577 |               95.60000 |                 0.21492 |     0.04798 |          0.98121 |               95.20000 |
| tmle      | glmnet_undersmooth_markov | Lasso           |       0.00197 |                 0 |                0.99224 |                     92.93333 | 0.00196 |           0 |          0.97520 |               95.60000 |                 0.21501 |     0.04807 |          0.98063 |               95.20000 |
| tmle      | glm                       | NA              |       0.00214 |                 0 |                1.05435 |                     91.86667 | 0.00203 |           0 |          0.96475 |               95.06667 |                 0.22094 |     0.04850 |          1.00320 |               95.60000 |
| tmle      | RF                        | NA              |       0.00319 |                 0 |                2.96372 |                     86.40000 | 0.00277 |           0 |          2.32506 |               92.00000 |                 0.28096 |     0.01065 |          2.72232 |               92.20000 |

Result: Ridge performs best

#### Comparison of penalized regression setup

All from Ridge models

<!-- = $\lambda$ selected at the minimum cross-validated SE -->
<!-- penalized regressions with $\lambda$ selected at the minimum cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L nodes in $g$ and $Q$ estimation. -->

| Undersmoothed        | Truncation  | Markov      | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:---------------------|:------------|:------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| Undersmoothed Lambda | g \< 0.01   | Markov L    |       0.00182 |             0e+00 |                0.82617 |                     93.46667 | 0.00184 |       0e+00 |          0.82992 |               95.20000 |                 0.20108 |     0.06054 |          0.81722 |               95.73333 |
| Undersmoothed Lambda | g \< 0.01   | All prior L |       0.00177 |             0e+00 |                0.80152 |                     93.86667 | 0.00186 |       0e+00 |          0.83601 |               94.66667 |                 0.20156 |     0.06268 |          0.80505 |               95.33333 |
| Undersmoothed Lambda | Untruncated | Markov L    |       0.00225 |             1e-05 |                0.68881 |                     93.37748 | 0.00227 |       1e-05 |          0.69410 |               95.89404 |                 0.24619 |     0.10833 |          0.74800 |               94.83444 |
| Undersmoothed Lambda | Untruncated | All prior L |       0.00227 |             1e-05 |                0.69077 |                     93.08511 | 0.00231 |       1e-05 |          0.70015 |               95.21277 |                 0.24826 |     0.10953 |          0.75014 |               94.81383 |
| Min SE Lambda        | Untruncated | Markov L    |       0.00195 |             0e+00 |                0.96458 |                     87.46667 | 0.00255 |       0e+00 |          1.20705 |               75.46667 |                 0.32056 |     0.08773 |          1.08226 |               83.06667 |
| Min SE Lambda        | Untruncated | All prior L |       0.00196 |             0e+00 |                0.94764 |                     88.13333 | 0.00258 |       0e+00 |          1.19545 |               76.00000 |                 0.32203 |     0.09015 |          1.07252 |               83.06667 |
| Min SE Lambda        | g \< 0.01   | Markov L    |       0.00197 |             0e+00 |                1.07245 |                     81.46667 | 0.00261 |       0e+00 |          1.35413 |               67.46667 |                 0.32959 |     0.07945 |          1.16926 |               78.66667 |
| Min SE Lambda        | g \< 0.01   | All prior L |       0.00197 |             0e+00 |                1.06886 |                     81.73333 | 0.00265 |       0e+00 |          1.36626 |               65.46667 |                 0.33201 |     0.08032 |          1.17147 |               78.40000 |

Result: Undersmoothed Ridge with default truncation and markov process
for L works best

## Comparison of IC and Bootstrap based coverages

    ## `summarise()` has grouped output by 'Target_parameter', 'sim_iter'. You can
    ## override using the `.groups` argument.

| Variance.estimator | Y\_.A.1..Coverage | RD.Coverage | RR.Coverage |
|:-------------------|------------------:|------------:|------------:|
| Influence curve    |          94.53333 |        90.0 |        94.4 |
| Bootstrap          |          94.20000 |        94.8 |        95.4 |

Results: Bootstrap has better coverage for RR and especially RD, but
slightly worse that for the mean of $Y_{A=1}$

## Plots of IC vs Bootstrap coverage

For the coverage of the risk difference (Red line is true RD)

<!-- # ```{r} -->
<!-- #  -->
<!-- # ic_plot_df <- res %>% filter(Estimator="tmle") -->
<!-- #  -->
<!-- # ggplot(ic_plot_df %>% -->
<!-- #          mutate(iteration=row_number()), aes(x=iteration)) + -->
<!-- #          geom_linerange(aes(ymin=boot.CI1, ymax=boot.CI2), alpha=0.5) + -->
<!-- #          coord_flip() + -->
<!-- #          geom_hline(yintercept = truth[10,4]) + ggtitle("IC") -->
<!-- #           -->
<!-- # ggplot(sim_res %>% filter(Target_parameter=="ATE") %>% arrange(boot.CI1) %>% ungroup() %>% -->
<!-- #          mutate(iteration=row_number()), aes(x=iteration)) + -->
<!-- #          geom_linerange(aes(ymin=boot.CI1, ymax=boot.CI2), alpha=0.5) + -->
<!-- #          coord_flip() + -->
<!-- #          geom_hline(yintercept = truth[10,4]) + ggtitle("Bootstrap") -->
<!-- #           -->
<!-- #         -->
<!-- #  -->
<!-- # ``` -->

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

## Updated needed in the registry analysis

Bochra was helping us rerun the analysis on the Danish registry based on
the updated simulation results with 500 iterations (undersmoothed ridge
won). Because of updated in (I think) TMLE_for_breakfast, there are now
D nodes for competing events, and when we run the old code we get the
error:

“Cannot both specify determistic.Q.nodes and Dnodes”

I think I just need to update both the run_tmle() and prepare_tmle()
functions to be able to pass a Dnodes=NULL argument through the
functions, because right now Dnodes are automatically set to be any
variable with “dead” within the name. But do I do that on Statistics
Denmark or make the changes on Ltmle_for_breakfast, and if the latter,
how do I export the update so it’s on Statistics Denmark?
