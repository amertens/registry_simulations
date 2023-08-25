Diabetes dementia analysis simulation results
================
Andrew Mertens

## Simulation setup

- N: 100,000

- Iterations: 500

- True $Y_{\bar{A}=1}$: `R truth[10,2]`

- True RD: `R truth[10,4]`

- True RR: `R truth[10,3]`

## Simulation results

    ##  [1] "X"                      "Estimator"              "estimator"             
    ##  [4] "N_reps"                 "abs_bias_Ya1"           "estimator_variance_Ya1"
    ##  [7] "mean_variance_Ya1"      "bias_se_ratio_Ya1"      "coverage_Ya1"          
    ## [10] "O_coverage_Ya1"         "abs_bias_RD"            "estimator_variance_RD" 
    ## [13] "mean_variance_RD"       "bias_se_ratio_RD"       "coverage_RD"           
    ## [16] "O_coverage_RD"          "abs_log_bias_RR"        "estimator_variance_RR" 
    ## [19] "mean_variance_RR"       "bias_se_ratio_RR"       "coverage_RR"           
    ## [22] "O_coverage_RR"

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

| Estimator | Algorithm | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:----------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| tmle      | EN        |       0.00182 |                 0 |                1.15564 |                         88.0 | 0.00243 |           0 |          1.45972 |                   74.6 |                 0.29780 |     0.05324 |          1.29068 |                   83.8 |
| tmle      | glm       |       0.00226 |                 0 |                1.09648 |                         91.8 | 0.00207 |           0 |          0.96961 |                   95.8 |                 0.22751 |     0.04808 |          1.03758 |                   95.4 |
| tmle      | glmnet    |       0.00181 |                 0 |                1.16584 |                         88.0 | 0.00244 |           0 |          1.48472 |                   74.6 |                 0.29728 |     0.05170 |          1.30745 |                   83.4 |
| tmle      | RF        |       0.00319 |                 0 |                2.96372 |                         86.4 | 0.00277 |           0 |          2.32506 |                   92.0 |                 0.28096 |     0.01065 |          2.72232 |                   92.2 |
| tmle      | ridge     |       0.00189 |                 0 |                1.00184 |                         84.4 | 0.00253 |           0 |          1.28295 |                   69.2 |                 0.31387 |     0.08017 |          1.10851 |                   80.2 |

Result: Ridge performs best

#### Comparison of penalized regression setup

Define different options

All from Lasso models

= $\lambda$ selected at the minimum cross-validated SE

penalized regressions with $\lambda$ selected at the minimum
cross-validated SE, truncation of $g$ at \< 0.01, and using all prior L
nodes in $g$ and $Q$ estimation.

| Estimator | Algorithm                            | Y\_{A=1} bias | Y\_{A=1} variance | Y\_{A=1} bias SE ratio | Y\_{A=1} oracle 95% coverage | RD bias | RD variance | RD bias SE ratio | RD oracle 95% coverage | RR log-transformed bias | RR variance | RR bias SE ratio | RR oracle 95% coverage |
|:----------|:-------------------------------------|--------------:|------------------:|-----------------------:|-----------------------------:|--------:|------------:|-----------------:|-----------------------:|------------------------:|------------:|-----------------:|-----------------------:|
| iptw      | ridge                                |       0.00324 |             0e+00 |                1.70973 |                     70.80000 | 0.00426 |       0e+00 |          2.14167 |               55.40000 |                 0.16829 |     0.02933 |          0.98258 |               93.40000 |
| iptw      | ridge_markov                         |       0.00324 |             0e+00 |                1.71849 |                     70.80000 | 0.00426 |       0e+00 |          2.15277 |               55.40000 |                 0.16829 |     0.02903 |          0.98769 |               93.40000 |
| iptw      | ridge_markov_untruncated             |       0.00324 |             0e+00 |                1.56001 |                     70.80000 | 0.00426 |       0e+00 |          1.97011 |               55.40000 |                 0.16829 |     0.03501 |          0.89944 |               93.40000 |
| iptw      | ridge_undersmooth                    |       0.00324 |             1e-05 |                1.42939 |                     70.80000 | 0.00426 |       1e-05 |          1.86101 |               55.40000 |                 0.16829 |     0.04117 |          0.82938 |               93.40000 |
| iptw      | ridge_undersmooth_markov             |       0.00324 |             1e-05 |                1.43083 |                     70.80000 | 0.00426 |       1e-05 |          1.86325 |               55.40000 |                 0.16829 |     0.04106 |          0.83052 |               93.40000 |
| iptw      | ridge_undersmooth_markov_untruncated |       0.00324 |             1e-05 |                0.96640 |                     70.94188 | 0.00426 |       1e-05 |          1.26567 |               55.31062 |                 0.16859 |     0.08796 |          0.56843 |               93.38677 |
| iptw      | ridge_undersmooth_untruncated        |       0.00323 |             1e-05 |                0.94901 |                     70.74830 | 0.00427 |       1e-05 |          1.24725 |               55.32880 |                 0.16717 |     0.09098 |          0.55423 |               93.42404 |
| iptw      | ridge_untruncated                    |       0.00324 |             0e+00 |                1.52670 |                     70.80000 | 0.00426 |       0e+00 |          1.93036 |               55.40000 |                 0.16829 |     0.03654 |          0.88037 |               93.40000 |
| tmle      | ridge                                |       0.00189 |             0e+00 |                1.00184 |                     84.40000 | 0.00253 |       0e+00 |          1.28295 |               69.20000 |                 0.31387 |     0.08017 |          1.10851 |               80.20000 |
| tmle      | ridge_markov                         |       0.00188 |             0e+00 |                1.00531 |                     84.20000 | 0.00249 |       0e+00 |          1.26925 |               70.60000 |                 0.31116 |     0.07926 |          1.10523 |               80.60000 |
| tmle      | ridge_markov_untruncated             |       0.00185 |             0e+00 |                0.89516 |                     88.60000 | 0.00240 |       0e+00 |          1.11740 |               77.00000 |                 0.29864 |     0.08773 |          1.00822 |               84.60000 |
| tmle      | ridge_undersmooth                    |       0.00186 |             1e-05 |                0.83031 |                     94.60000 | 0.00185 |       1e-05 |          0.81642 |               96.20000 |                 0.20352 |     0.06227 |          0.81556 |               95.20000 |
| tmle      | ridge_undersmooth_markov             |       0.00194 |             1e-05 |                0.86738 |                     93.60000 | 0.00186 |       1e-05 |          0.82423 |               96.20000 |                 0.20619 |     0.05989 |          0.84253 |               96.00000 |
| tmle      | ridge_undersmooth_markov_untruncated |       0.00236 |             1e-05 |                0.71299 |                     92.98597 | 0.00226 |       1e-05 |          0.68211 |               96.19238 |                 0.24777 |     0.10804 |          0.75380 |               94.78958 |
| tmle      | ridge_undersmooth_untruncated        |       0.00231 |             1e-05 |                0.68792 |                     93.19728 | 0.00222 |       1e-05 |          0.66002 |               95.91837 |                 0.23960 |     0.11118 |          0.71859 |               95.46485 |
| tmle      | ridge_untruncated                    |       0.00185 |             0e+00 |                0.87856 |                     89.60000 | 0.00242 |       0e+00 |          1.10716 |               77.40000 |                 0.30002 |     0.09028 |          0.99853 |               84.80000 |

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

<!-- # ```{r,  echo=F} -->
<!-- # # for(x in 1:nrow(d)){ -->
<!-- # #   #they just need to follow until they get the outcome or until censored (take min) -->
<!-- # #   max_fu <- min(match(1, data[x,..Ynodes]),match(1, data[x,..Cnodes]),N_time,na.rm=T) -->
<!-- # #   a_sub <- Anodes[1:max_fu] -->
<!-- # #   data$regime_fu <- max_fu -->
<!-- # #   #if sum of subsetted A nodes is = to follw up time, then they followed full regime -->
<!-- # #   data$full_regime_1[[x]] <- ifelse(sum(data[x,..a_sub])==max_fu,TRUE, FALSE) -->
<!-- # #   #if sum of subsetted a nodes is 0 then they had no glp1 throughout follow up -->
<!-- # #   data$full_regime_0[[x]] <- ifelse(sum(data[x,..a_sub])==0,TRUE, FALSE) -->
<!-- # # } -->
<!-- #  -->
<!-- #  -->
<!-- #   data[,maxfu_y:= N_time-rowSums(.SD),.SDcols=c(Ynodes)]#max fu y -->
<!-- #   data[,maxfu_c:= N_time-rowSums(.SD),.SDcols=c(Cnodes)]#max fu c -->
<!-- #   data[,maxfu_total:= min(N_time,maxfu_y,maxfu_c)] #total time fu (taking min) -->
<!-- #   a_sub <- Anodes[1:max_fu] -->
<!-- #   data$regime_fu <- max_fu -->
<!-- #   data[,a_sub_sum:= rowSums(.SD),.SDcols=a_sub]#manodes follow in fu -->
<!-- #  -->
<!-- #   #if sum of subsetted A nodes is = to follw up time, then they followed full regime -->
<!-- #   data[a_sub_sum==max_fu,full_regime_1:=T,] -->
<!-- #    -->
<!-- #   #if sum of subsetted a nodes is 0 then they had no glp1 throughout follow up -->
<!-- #   data$full_regime_0[[x]] <- ifelse(sum(data[x,..a_sub])==0,TRUE, FALSE) -->
<!-- # } -->
<!-- # data$switchers<- ifelse(data$full_regime_1==F & data$full_regime_0==F, TRUE,FALSE) -->
<!-- #  -->
<!-- # tab <- CreateCatTable(c("full_regime_1","full_regime_0","switchers","regime_fu"),data=data) -->
<!-- # kableone(tab) -->
<!-- #  -->
<!-- #  -->
<!-- # ``` -->
<!-- #  -->
<!-- #  -->
<!-- # ## Outcome -->
<!-- #  -->
<!-- # ```{r} -->
<!-- #  -->
<!-- # tab <- CreateCatTable(c(Ynodes),data=data) -->
<!-- # kableone(tab) -->
<!-- #  -->
<!-- #  -->
<!-- # ``` -->
<!-- #  -->
<!-- #  -->
<!-- #  -->
<!-- # ## Regime support stratified by outcome (overall) -->
<!-- #  -->
<!-- # ```{r} -->
<!-- #  -->
<!-- # tab <- CreateCatTable(c("full_regime_1","full_regime_0","switchers"),data=data, strata="event_dementia_10") -->
<!-- # kableone(tab) -->
<!-- #  -->
<!-- # ``` -->
<!-- #  -->
<!-- # ## Censoring -->
<!-- #  -->
<!-- # ```{r, echo=F} -->
<!-- #  -->
<!-- # # for(i in 1:(N_time-2)){ -->
<!-- # #   cens.node <- (paste0("censor_",(i+1))) -->
<!-- # #   data[get(paste0("censor_",i))==1,get(cens.node) := replace(.SD, .SD == 0, 1), .SDcols = cens.node] -->
<!-- # # } -->
<!-- # #  -->
<!-- # # d[sum_death < sum_dementia, (death.nodes) := replace(.SD, .SD == 1, 0), .SDcols = death.nodes] -->
<!-- # tab <- CreateCatTable(c(Cnodes),data=data) -->
<!-- # kableone(tab) -->
<!-- #  -->
<!-- # ``` -->
<!-- #  -->
<!-- #  -->
<!-- # ## Covariate characteristics by those following treatment rules at distinct time points -->
<!-- #  -->
<!-- # ```{r,include=F,echo=FALSE} -->
<!-- #  -->
<!-- # find_followup <- function(data,t){ -->
<!-- #     #subset to risk set (no outcome or censor at time t) -->
<!-- #     data2 <- data[get(paste0("censor_",t))!=1 & get(paste0("event_dementia_",t))!=1]  -->
<!-- #      -->
<!-- #   for(x in 1:nrow(data2)){ -->
<!-- #     #subset Anodes to check   -->
<!-- #     a_sub <- Anodes[1:t] -->
<!-- #     #if sum of subsetted A nodes is = to follw up time, then they followed full regime -->
<!-- #     data2$full_regime_1[[x]] <- ifelse(sum(data2[x,..a_sub])==t,TRUE, FALSE) -->
<!-- #     #if sum of subsetted a nodes is 0 then they had no glp1 throughout follow up -->
<!-- #     data2$full_regime_0[[x]] <- ifelse(sum(data2[x,..a_sub])==0,TRUE, FALSE) -->
<!-- #      -->
<!-- #   } -->
<!-- #    -->
<!-- #   names <- c("insulin_","any.malignancy_", "chronic.pulmonary.disease_", -->
<!-- #              "hypertension_", "myocardial.infarction_","ischemic.heart.disease_", "heart.failure_","renal.disease_","sglt2_inhib_", "event_dementia_") -->
<!-- #   cov_timedep <- outer(names,0:t,paste0) -->
<!-- #   #updated exposure var, drop folks who do not follow either rule for the purposes of the table -->
<!-- #   data_sub <- data2[data2$full_regime_1==T | data2$full_regime_0==T,] -->
<!-- #   data_sub$exposure_regime <- ifelse(data_sub$full_regime_1==T,"Full GLP1", "No GLP1") -->
<!-- #   cov_set <- c("sex",cov_timedep) -->
<!-- #  -->
<!-- #  -->
<!-- # return(list(cov_set,data_sub)) -->
<!-- # } -->
<!-- #  -->
<!-- #  -->
<!-- # # kableone(find_followup(data=data,t=2)) -->
<!-- #  -->
<!-- # ``` -->
<!-- #  -->
<!-- #  -->
<!-- # ####  t=2 -->
<!-- #  -->
<!-- # ```{r,echo=FALSE} -->
<!-- #  -->
<!-- # t <- find_followup(data=data,t=2) -->
<!-- # tab <- CreateCatTable(t[[1]],data=t[[2]], strata="exposure_regime") -->
<!-- # kableone(tab) -->
<!-- # ``` -->
<!-- #  -->
<!-- # #### t=3 -->
<!-- #  -->
<!-- # ```{r,echo=FALSE} -->
<!-- #  -->
<!-- # t <- find_followup(data=data,t=3) -->
<!-- # tab <- CreateCatTable(t[[1]],data=t[[2]], strata="exposure_regime") -->
<!-- # kableone(tab) -->
<!-- # ``` -->
<!-- #  -->
<!-- # ## Unbounded g weights from LTMLE run, TCP 1 -->
<!-- #  -->
<!-- #  -->
<!-- # ```{r, echo=FALSE,warning=F} -->
<!-- # load("../data/NOTRANSFER_glp1_any_static11.RData") -->
<!-- # cum.g.unbounded <- res_RR$cum.g.unbounded -->
<!-- #  -->
<!-- # tx1 <- cbind(cum.g.unbounded[,,1],rep(1,nrow(cum.g.unbounded[,,1]))) -->
<!-- # tx0 <- cbind(cum.g.unbounded[,,2],rep(0,nrow(cum.g.unbounded[,,2]))) -->
<!-- #  -->
<!-- # all<-as.data.frame(rbind(tx1,tx0)) -->
<!-- #  -->
<!-- # names(all)<- c("A1","C1","A2","C2","A3","C3","A4","C4", "A5","C5", -->
<!-- #                "A6","C6","A7","C7","A8","C8","A9","C9","A10","C10","regime") -->
<!-- # gathered <- tidyr::gather(data.frame(all),key="node",value="value",1:20) -->
<!-- # gathered$regime <- as.character(gathered$regime) -->
<!-- #  -->
<!-- # gather_sub <- gathered[gathered$node%in%c(names(all)[grep("A",names(all))]),] -->
<!-- #  -->
<!-- # library(ggplot2) -->
<!-- # ggplot(gather_sub, aes(x=value, fill=regime)) + -->
<!-- #   geom_histogram(position="identity", alpha=0.2)+  -->
<!-- #   theme_classic() + -->
<!-- #   facet_wrap((~node))+theme(legend.position ="right") -->
<!-- # ``` -->
