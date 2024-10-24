

rm(list=ls())
load(paste0(here::here(),"/data/sim_results/sim_res_markov.Rdata"))

gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)
library(targets)
library(tarchetypes)
library(parallel)

lapply(c("fst","lava","ltmle","data.table","tidyverse","glmnet","Matrix","matrixStats","speedglm","parallel","caret","foreach","clustermq"), FUN = function(X) {
  do.call("require", list(X)) 
})

# -------------------------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------------------------

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
set.seed(12345)

# # Set the simulation hyperparameters
#simulated dataset size
n_df=100000 
#set time horizon
time=10
#set
#Longitudinal variables possibly following the markov process
Markov_variables=c("heart.failure","renal.disease","chronic.pulmonary.disease", "any.malignancy"  ,         
                   "ischemic.heart.disease","myocardial.infarction","hypertension","stroke" ,                  
                   "bb","ccb","rasi","thiazid",
                   "loop","mra","copd_med"  )

df = readRDS(paste0(here::here(),"/data/sim_data/simulated_data_",1,".rds"))


  library="glm"
           SL.Control=NULL
           dataset_num=1
           n=1
           time=2
           gbounds=c(0.1,1)
           null_sim=FALSE
           n_bootstrap_samples=0
           Markov_variables=Markov_variables
           tmle_var=FALSE

    nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
    nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)
    
    # model <- targets::tar_read_raw("lava_model")
    # simulated_data = simulate_data(lava_model = model, n = n)
    # simulated_data = clean_sim_data(simulated_data, N_time=time)
    simulated_data = readRDS(paste0(here::here(),"/data/sim_data/simulated_data_",dataset_num,".rds"))
    simulated_data = simulated_data[1:10000,]
    

    # # #set up analysis:
    simulated_data_list <- get_simulated_data_list(simulated_data = simulated_data, time_horizon=time)
    #simulated_data_list$outcome_data
    # res = run_ltmle(name_outcome="dementia",
    #                 time_horizon=time,
    #                 test=FALSE,
    #                 outcome_data=simulated_data_list$outcome_data,
    #                 regimen_data=simulated_data_list$regimen_data,
    #                 baseline_data=simulated_data_list$sim_baseline_covariates,
    #                 timevar_data=simulated_data_list$sim_time_covariates,
    #                 det.Q.function=NULL,# now build-in
    #                 gbounds=gbounds,
    #                 B_bootstrap_samples=n_bootstrap_samples,
    #                 SL.library=library,
    #                 Markov=Markov_variables,
    #                 SL.cvControl=SL.Control,
    #                 tmle_var=tmle_var,
    #                 verbose=TRUE)
    
    
    # res=res[[1]][[1]]
    # res=res[1]
    # #res
    # 
    # output=summary(res$Ltmle_fit)
    # output_iptw= summary(res$Ltmle_fit, estimator="iptw")
    # output=bind_rows(output, output_iptw)
    # output$dataset_num <- dataset_num
    # output
    
    name_outcome="dementia"
                    time_horizon=time
                    test=FALSE
                    outcome_data=simulated_data_list$outcome_data
                    regimen_data=simulated_data_list$regimen_data
                    baseline_data=simulated_data_list$sim_baseline_covariates
                    timevar_data=simulated_data_list$sim_time_covariates
                    det.Q.function=NULL
                    B_bootstrap_samples=n_bootstrap_samples
                    SL.library=library
                    Markov=Markov_variables
                    SL.cvControl=SL.Control
                    verbose=TRUE
    

             sub_set=NULL
             test=FALSE
             gbounds=c(0.01,1)
           
             # SL.library="glmnet"
             # SL.cvControl=list(selector="undersmoothed",alpha=1)
             B_bootstrap_samples=0
             seeds=NULL
             verbose=FALSE
             reduce=TRUE
             tmle_var=FALSE
      require(foreach,quietly=TRUE)
      require(data.table,quietly=TRUE)
             
        tk=2
      #result <- foreach(tk=time_horizon)%do%{
        abar <- list(rep(1,tk),rep(0,tk))
        REG = names(regimen_data)[1]
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
          
          #add truncation
          pl$gbounds = gbounds
          if(tmle_var){
            pl$deterministic.Q.function=deterministic.Q.function=NULL
            pl$variance.method = "tmle"
            pl$data[is.na(pl$data)] <- 0
            pl$data[is.na(pl$data)] <- "censored"
          }
          
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
          
          pl$iptw.only<-FALSE
          tryfit <- try(fit <- do.call(Ltmle,pl))
          
          save.image(here::here("my_environment.RData"))
          
          
          pl4 <- pl2 <- pl
          pl2$SL.library <- "glmnet"
          pl2$SL.cvControl <- list(selector="undersmoothed",alpha=1)
          tryfit2 <- try(fit <- do.call(Ltmle,pl2))
          summary(tryfit, estimator="iptw")
          
          tryfit$fit$g[[1]]$GLP1RA_1
          tryfit2$fit$g[[1]]$GLP1RA_1
          
          summary(tryfit, estimator="iptw")
          summary(tryfit2, estimator="iptw")
          
          summary(tryfit)
          summary(tryfit2)
          
          pl3 <- pl
          pl3$info<-NULL
          tryfit3 <- try(fit <- do.call(ltmle,pl3))
          summary(tryfit3, estimator="iptw")

          pl4$iptw.only<-TRUE
          tryfit4 <- try(fit <- do.call(Ltmle,pl4))
          summary(tryfit4, estimator="iptw")
          
          pl5 <- pl
          pl5$time_horizon <- 2
          tryfit5 <- try(fit <- do.call(Ltmle2,pl5))
          
          
          

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
      #   names(loop)=names(regimen_data)[1:length(loop)]
      #   loop
      # }
      names(result)=paste0("time_horizon_",time_horizon)
      result
    #}