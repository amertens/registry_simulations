Diabetes dementia analysis simulation results
================
Andrew Mertens

## Simulation setup

- N: 100,000

- Iterations: 500

- True $Y_{\bar{A}=1}$: 0.0082798

- True RD: -0.0071862

- True RR: 0.5356286

## Simulation results

#### Comparison of IPTW and TMLE estimators

Both from Lasso models with $\lambda$ selected at the minimum
cross-validated SE

| Estimator | Algorithm | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| iptw      | glmnet    |       0.00324 |                 0 |                2.08002 |                         70.8 | 0.00426 |           0 |          2.57546 |                   55.4 |                 0.16829 |     0.01981 |          1.19573 |                   93.4 |
| tmle      | glmnet    |       0.00181 |                 0 |                1.16584 |                         88.0 | 0.00244 |           0 |          1.48472 |                   74.6 |                 0.29728 |     0.05170 |          1.30745 |                   83.4 |

Result: TMLE performs better

#### Comparison of algorithms

All penalized regressions with $\lambda$ selected at the minimum
cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L
nodes in $g$ and $Q$ estimation.

RF= random forest

| Estimator | Algorithm Alpha | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | Elastic Net     |       0.00182 |                 0 |                1.15564 |                         88.0 | 0.00243 |           0 |          1.45972 |                   74.6 |                 0.29780 |     0.05324 |          1.29068 |                   83.8 |
| tmle      | NA              |       0.00226 |                 0 |                1.09648 |                         91.8 | 0.00207 |           0 |          0.96961 |                   95.8 |                 0.22751 |     0.04808 |          1.03758 |                   95.4 |
| tmle      | Lasso           |       0.00181 |                 0 |                1.16584 |                         88.0 | 0.00244 |           0 |          1.48472 |                   74.6 |                 0.29728 |     0.05170 |          1.30745 |                   83.4 |
| tmle      | NA              |       0.00319 |                 0 |                2.96372 |                         86.4 | 0.00277 |           0 |          2.32506 |                   92.0 |                 0.28096 |     0.01065 |          2.72232 |                   92.2 |
| tmle      | Ridge           |       0.00189 |                 0 |                1.00184 |                         84.4 | 0.00253 |           0 |          1.28295 |                   69.2 |                 0.31387 |     0.08017 |          1.10851 |                   80.2 |

Result: Ridge performs best

#### Comparison of penalized regression setup

All from Ridge models

<!-- = $\lambda$ selected at the minimum cross-validated SE -->
<!-- penalized regressions with $\lambda$ selected at the minimum cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L nodes in $g$ and $Q$ estimation. -->

| Estimator | Truncation  | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | g \< 0.01   |       0.00189 |             0e+00 |                1.00184 |                     84.40000 | 0.00253 |       0e+00 |          1.28295 |               69.20000 |                 0.31387 |     0.08017 |          1.10851 |               80.20000 |
| tmle      | g \< 0.01   |       0.00188 |             0e+00 |                1.00531 |                     84.20000 | 0.00249 |       0e+00 |          1.26925 |               70.60000 |                 0.31116 |     0.07926 |          1.10523 |               80.60000 |
| tmle      | Untruncated |       0.00185 |             0e+00 |                0.89516 |                     88.60000 | 0.00240 |       0e+00 |          1.11740 |               77.00000 |                 0.29864 |     0.08773 |          1.00822 |               84.60000 |
| tmle      | g \< 0.01   |       0.00186 |             1e-05 |                0.83031 |                     94.60000 | 0.00185 |       1e-05 |          0.81642 |               96.20000 |                 0.20352 |     0.06227 |          0.81556 |               95.20000 |
| tmle      | g \< 0.01   |       0.00194 |             1e-05 |                0.86738 |                     93.60000 | 0.00186 |       1e-05 |          0.82423 |               96.20000 |                 0.20619 |     0.05989 |          0.84253 |               96.00000 |
| tmle      | Untruncated |       0.00236 |             1e-05 |                0.71299 |                     92.98597 | 0.00226 |       1e-05 |          0.68211 |               96.19238 |                 0.24777 |     0.10804 |          0.75380 |               94.78958 |
| tmle      | Untruncated |       0.00231 |             1e-05 |                0.68792 |                     93.19728 | 0.00222 |       1e-05 |          0.66002 |               95.91837 |                 0.23960 |     0.11118 |          0.71859 |               95.46485 |
| tmle      | Untruncated |       0.00185 |             0e+00 |                0.87856 |                     89.60000 | 0.00242 |       0e+00 |          1.10716 |               77.40000 |                 0.30002 |     0.09028 |          0.99853 |               84.80000 |

Result: Undersmoothed Ridge with default truncation and markov process
for L works best

## Comparison of IC and Bootstrap based coverages

    ## `summarise()` has grouped output by 'Target_parameter', 'sim_iter'. You can
    ## override using the `.groups` argument.

| Variance.estimator | Y\_.A.1..Coverage | RD.Coverage | RR.Coverage |
|:-------------------|------------------:|------------:|------------:|
| Influence curve    |              95.2 |        91.4 |        94.2 |
| Bootstrap          |              94.2 |        94.8 |        95.4 |

Results: Bootstrap has better coverage for RR and especially RD, but
slightly worse that for the mean of $Y_{A=1}$

## Boxplots

- Make two boxplots of bootstrap vs IC

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
