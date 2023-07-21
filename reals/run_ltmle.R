run_ltmle <- function(name_outcome,
                      time_horizon,
                      sub_set=NULL,
                      test=FALSE,
                      regimen_data,
                      outcome_data,
                      baseline_data,
                      timevar_data,
                      Markov=NULL,
                      det.Q.function,
                      SL.library="glmnet",
                      SL.cvControl=list(selector="undersmoothed",alpha=1),
                      B_bootstrap_samples=0,
                      seeds,
                      verbose=FALSE,
                      reduce=TRUE){
    require(foreach,quietly=TRUE)
    require(data.table,quietly=TRUE)
    result <- foreach(tk=time_horizon)%do%{
        abar <- list(rep(1,tk),rep(0,tk))
        loop <- foreach(REG = names(regimen_data))%do%{
            bsl_covariates <- baseline_data
            setkey(bsl_covariates,pnr)
            ## add baseline adjustment to subset analysis
            if (length(sub_set)>0 & length(sub_set$adj)>0){
                sdat=sub_set$data[,c("pnr",sub_set$adj),with=FALSE]
                setkey(sdat,pnr)
                bsl_covariates <- sdat[bsl_covariates]
            }
            if (length(sub_set)>0){
                sub_id <- sub_set$data[["pnr"]]
                if(length(sub_id)==0)stop("No data in subset defined by variable: ",sub_set$var)
            } else{
                sub_id <- NULL
            }
            regimens <- REG
            # all timevarying covariates but not treatment
            # should only enter the formula with their last value
            if (length(Markov)>0){
                markov <- sub("_0","",grep("_0",names(timevar_data),value=TRUE))
                if (is.character(Markov)){
                    markov = intersect(markov,Markov)
                }
            } else{
                markov=""
            }
            pl=prepare_Ltmle(regimen_data=regimen_data[[REG]],
                             outcome_data=outcome_data,
                             name_outcome=name_outcome,
                             name_regimen=regimens,
                             name_censoring = "Censored",
                             censored_label = "censored",
                             name_comp.event = "Dead",
                             baseline_data=bsl_covariates,
                             timevar_data=timevar_data,
                             time_horizon=tk,
                             subset_id=sub_id,
                             test=test,
                             SL.library=SL.library,
                             Markov=markov,
                             deterministic.Q.function=det.Q.function,
                             abar=abar)
            if (verbose){
                cat("Run Ltmle for regimen ",
                    paste0(regimens,collapse=","),
                    " and outcome ",
                    name_outcome,
                    "\n",
                    sep="")
                pl$verbose <- TRUE}
            if (length(SL.cvControl)>0)
                pl$SL.cvControl <- SL.cvControl
            if (B_bootstrap_samples>0){
                if(missing(seeds))
                    seeds <- sample(1:1000000,size=B_bootstrap_samples)
                library(doParallel)
                library(foreach)
                library(data.table)
                message("Bootstrapping")
                tb=txtProgressBar(max=B_bootstrap_samples,width=20,style=3)
                bootfit <- foreach(b=1:B_bootstrap_samples,
                                   .combine="rbind",.export=c("run_ltmle","summary.Ltmle","prepare_Ltmle","Ltmle","get_formulas","get_ltmle_data","get_rhs","get_subset_data","merge_data"),
                                   .packages=c("data.table"))%dopar%{
                                       setTxtProgressBar(pb=tb,b)
                                       pl.b <- data.table::copy(pl)
                                       pl.b$verbose=FALSE
                                       set.seed(seeds[[b]])
                                       pl.b$data <- pl$data[sample(1:.N,replace=TRUE,size=.N)]
                                       pl.b$id=pl.b$data$pnr
                                       tryfit <- try(fit.b <- do.call(Ltmle,pl.b))
                                       if (inherits(tryfit,"try-error")) return(NULL)
                                       est.b=summary(fit.b)[,.(Target_parameter=Target_parameter,estimate=estimate,b=b)]
                                       est.b
                                   }
            }
            if (verbose)print(paste0("Fitting Ltmle"," ",REG))
            tryfit <- try(fit <- do.call(Ltmle,pl))
            if (inherits(tryfit,"try-error")) return(NULL) # browser()
            if (reduce){
                fit$call <- NULL
                fit$cum.g <- fit$cum.g.used <- fit$cum.g.unbounded <- NULL
                ## fit$IC <- NULL
                fit$Qstar <- NULL
            }
            x=c(list(Ltmle_fit=fit,time_horizon=tk,regimen=REG),
                # formula are potential data/environment collectors
                # when object is saved hence we not include them
                # in the output
                ## Qform=Qform,
                ## gform=gform,
                with(pl,list(Anodes=Anodes,
                             Cnodes=Cnodes,
                             Lnodes=Lnodes,
                             Ynodes=Ynodes,
                             abar=abar,
                             deterministic.Q.function=deterministic.Q.function,
                             SL.library=SL.library,
                             SL.cvControl=SL.cvControl)))
            if (B_bootstrap_samples>0)
                x$bootfit <- bootfit
            x
        }
        names(loop)=names(regimen_data)[1:length(loop)]
        loop
    }
    names(result)=paste0("time_horizon_",time_horizon)
    result
}
