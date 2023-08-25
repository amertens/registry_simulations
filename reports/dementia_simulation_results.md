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

    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.
    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.
    ## `summarise()` has grouped output by 'Estimator'. You can override using the
    ## `.groups` argument.

## Simulation results

#### Comparison of IPTW and TMLE estimators

Both from Lasso models with $\lambda$ selected at the minimum
cross-validated SE. All tables below are sorted by the RD bias.

| Estimator | Algorithm | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | glmnet    |       0.00184 |                 0 |                1.19893 |                         86.9 | 0.00249 |           0 |          1.52618 |                   74.0 |                 0.30534 |     0.05190 |          1.34024 |                   83.1 |
| iptw      | glmnet    |       0.00324 |                 0 |                2.08002 |                         70.8 | 0.00426 |           0 |          2.57546 |                   55.4 |                 0.16829 |     0.01981 |          1.19573 |                   93.4 |

Result: TMLE performs better

#### Comparison of algorithms

All penalized regressions with $\lambda$ selected with undersmoothing,
truncation of $g$ at \< 0.01, and using the t-1 prior L nodes in $g$ and
$Q$ estimation.

RF= random forest

| Estimator | Algorithm                 | Algorithm Alpha | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:--------------------------|:----------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | ridge_undersmooth_markov  | Ridge           |       0.00188 |                 0 |                0.84774 |                         93.3 | 0.00186 |       1e-05 |          0.82992 |                   96.1 |                 0.20412 |     0.06025 |          0.83158 |                   96.0 |
| tmle      | EN_undersmooth_markov     | Elastic Net     |       0.00203 |                 0 |                1.01708 |                         92.7 | 0.00198 |       0e+00 |          0.97798 |                   95.9 |                 0.21796 |     0.04768 |          0.99822 |                   95.6 |
| tmle      | glmnet_undersmooth_markov | Lasso           |       0.00203 |                 0 |                1.01669 |                         92.9 | 0.00198 |       0e+00 |          0.97747 |                   95.9 |                 0.21807 |     0.04777 |          0.99768 |                   95.6 |
| tmle      | glm                       | NA              |       0.00220 |                 0 |                1.07738 |                         92.0 | 0.00206 |       0e+00 |          0.97150 |                   95.4 |                 0.22486 |     0.04822 |          1.02396 |                   95.7 |
| tmle      | RF                        | NA              |       0.00319 |                 0 |                2.96372 |                         86.4 | 0.00277 |       0e+00 |          2.32506 |                   92.0 |                 0.28096 |     0.01065 |          2.72232 |                   92.2 |

Result: Ridge performs best

#### Comparison of penalized regression setup

All from Ridge models

<!-- = $\lambda$ selected at the minimum cross-validated SE -->
<!-- penalized regressions with $\lambda$ selected at the minimum cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L nodes in $g$ and $Q$ estimation. -->

| Undersmoothed        | Truncation  | Markov      | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:---------------------|:------------|:------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| Undersmoothed Lambda | g \< 0.01   | Markov L    |       0.00188 |             0e+00 |                0.84774 |                     93.30000 | 0.00186 |       1e-05 |          0.82992 |               96.10000 |                 0.20412 |     0.06025 |          0.83158 |               96.00000 |
| Undersmoothed Lambda | g \< 0.01   | All prior L |       0.00182 |             0e+00 |                0.81819 |                     94.20000 | 0.00186 |       1e-05 |          0.82993 |               95.20000 |                 0.20325 |     0.06245 |          0.81335 |               95.40000 |
| Undersmoothed Lambda | Untruncated | Markov L    |       0.00229 |             1e-05 |                0.69915 |                     93.35317 | 0.00226 |       1e-05 |          0.68722 |               96.13095 |                 0.24666 |     0.10853 |          0.74872 |               94.94048 |
| Undersmoothed Lambda | Untruncated | All prior L |       0.00229 |             1e-05 |                0.69113 |                     93.22034 | 0.00228 |       1e-05 |          0.68574 |               95.55085 |                 0.24603 |     0.11046 |          0.74028 |               95.02119 |
| Min SE Lambda        | Untruncated | Markov L    |       0.00190 |             0e+00 |                0.92971 |                     87.90000 | 0.00247 |       0e+00 |          1.16188 |               76.70000 |                 0.31022 |     0.08788 |          1.04643 |               84.00000 |
| Min SE Lambda        | Untruncated | All prior L |       0.00191 |             0e+00 |                0.91316 |                     88.80000 | 0.00250 |       0e+00 |          1.15055 |               77.20000 |                 0.31154 |     0.09037 |          1.03632 |               84.00000 |
| Min SE Lambda        | g \< 0.01   | Markov L    |       0.00193 |             0e+00 |                1.03908 |                     82.80000 | 0.00255 |       0e+00 |          1.31246 |               69.60000 |                 0.32078 |     0.07935 |          1.13874 |               80.30000 |
| Min SE Lambda        | g \< 0.01   | All prior L |       0.00193 |             0e+00 |                1.03628 |                     83.20000 | 0.00259 |       0e+00 |          1.32561 |               67.70000 |                 0.32335 |     0.08021 |          1.14175 |               80.10000 |

Result: Undersmoothed Ridge with default truncation and markov process
for L works best

## Comparison of IC and Bootstrap based coverages

    ## `summarise()` has grouped output by 'Target_parameter', 'sim_iter'. You can
    ## override using the `.groups` argument.

| Variance.estimator | Y\_.A.1..Coverage | RD.Coverage | RR.Coverage |
|:-------------------|------------------:|------------:|------------:|
| Influence curve    |              94.7 |        90.3 |        94.2 |
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
