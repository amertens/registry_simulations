

get_simulated_data_list <- function(simulated_data, time_horizon=10){
  simulated_data[,pnr:=1:.N]
  names_time_covariates = c("heart.failure","renal.disease","chronic.pulmonary.disease","any.malignancy","ischemic.heart.disease","myocardial.infarction","hypertension","stroke","bb","ccb","rasi","thiazid","loop","mra","copd_med")
  sim_time_covariates <- simulated_data[,c("pnr",c(sapply(names_time_covariates,paste0,"_",0:(time_horizon-1)))),with = FALSE]
  sim_outcome<- simulated_data[,grep("pnr|dementia_|Censored|Dead", names(simulated_data)), with = FALSE]
  ## fix outcome to value 1 after first occurrence
  for (i in 1:time_horizon){
    for (j in ((i+1):time_horizon)){
      set(sim_outcome,i=which(sim_outcome[[paste0("dementia_",i)]]==1),j=paste0("dementia_",j),value=1)
    }
  }
  sim_baseline_covariates=simulated_data[,c("pnr","sex","agegroups","education","income","diabetes_duration"),with=FALSE]
  list(regimen_data = list("GLP1RA" = simulated_data[,grep("pnr|GLP1RA", names(simulated_data)), with = FALSE]),
       outcome_data = list(dementia=sim_outcome),
       sim_baseline_covariates = sim_baseline_covariates,
       sim_time_covariates = sim_time_covariates)
}









run_ltmle <- function(name_outcome,
                      time_horizon,
                      sub_set=NULL,
                      test=FALSE,
                      regimen_data,
                      outcome_data,
                      baseline_data,
                      timevar_data,
                      gbounds=c(0.01,1),
                      Markov=NULL,
                      det.Q.function,
                      SL.library="glmnet",
                      SL.cvControl=list(selector="undersmoothed",alpha=1),
                      B_bootstrap_samples=0,
                      seeds,
                      verbose=FALSE,
                      reduce=TRUE,
                      tmle_var=FALSE){
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

AsMatrix <-
  function (x) 
  {
    if (is.matrix(x)) {
      return(x)
    }
    else if (is.vector(x)) {
      dim(x) <- c(length(x), 1)
      return(x)
    }
    else {
      stop("AsMatrix input should be a matrix or vector")
    }
  }

Bound <-
  function (x, bounds) 
  {
    stopifnot(length(bounds) == 2 && !anyNA(bounds))
    x[x < min(bounds)] <- min(bounds)
    x[x > max(bounds)] <- max(bounds)
    return(x)
  }
CalcCumG <-
  function (g, gbounds) 
  {
    cum.g <- rowCumprods(g)
    return(list(unbounded = cum.g, bounded = Bound(cum.g, gbounds)))
  }
CalcG <-
  function (prob.A.is.1, cur.abar, deterministic.newdata) 
  {
    g <- matrix(NA_real_, nrow(prob.A.is.1), ncol(prob.A.is.1))
    g[!is.na(cur.abar) & cur.abar == 1] <- prob.A.is.1[!is.na(cur.abar) & 
                                                         cur.abar == 1]
    g[!is.na(cur.abar) & cur.abar == 0] <- 1 - prob.A.is.1[!is.na(cur.abar) & 
                                                             cur.abar == 0]
    g[deterministic.newdata] <- 1
    return(g)
  }
CalcGUnboundedToBoundedRatio <-
  function (g.list, nodes, final.Ynodes) 
  {
    CalcForFinalYNode <- function(num.AC.nodes) {
      if (num.AC.nodes == 0) 
        return(1)
      if (!anyNA(g.list$cum.g)) 
        return(AsMatrix(g.list$cum.g.unbounded[, num.AC.nodes, 
        ]/g.list$cum.g[, num.AC.nodes, ]))
      g.ratio1 <- matrix(NA, n, num.regimes)
      for (i in 1:num.regimes) {
        g.ratio.temp <- cbind(g.list$cum.g.meanL.unbounded[, 
                                                           num.AC.nodes, i, ]/g.list$cum.g.meanL[, num.AC.nodes, 
                                                                                                 i, ], g.list$cum.g.unbounded[, num.AC.nodes, 
                                                                                                                              i]/g.list$cum.g[, num.AC.nodes, i])
        index <- max.col(!is.na(g.ratio.temp), "last")
        g.ratio1[, i] <- g.ratio.temp[sub2ind(1:n, col = index, 
                                              num.rows = n)]
      }
      return(g.ratio1)
    }
    n <- dim(g.list$cum.g)[1]
    num.regimes <- dim(g.list$cum.g)[3]
    num.final.Ynodes <- length(final.Ynodes)
    g.ratio <- array(dim = c(n, num.regimes, num.final.Ynodes))
    for (j in 1:num.final.Ynodes) {
      num.AC.nodes <- sum(nodes$AC < final.Ynodes[j])
      g.ratio[, , j] <- CalcForFinalYNode(num.AC.nodes)
    }
    return(g.ratio)
  }
CalcIC <-
  function (Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, 
            regimes.with.positive.weight) 
  {
    n <- nrow(Qstar)
    num.regimes <- ncol(Qstar)
    num.betas <- dim(h.g.ratio)[3]
    IC <- matrix(0, nrow = n, ncol = num.betas)
    for (i in regimes.with.positive.weight) {
      index <- uncensored & intervention.match[, i]
      if (any(h.g.ratio[index, i, ] != 0)) {
        IC[index, ] <- IC[index, ] + (Qstar.kplus1[index, 
                                                   i] - Qstar[index, i]) * h.g.ratio[index, i, ]
      }
    }
    return(IC)
  }
CalcInterventionMatchArray <-
  function (data, regimes, Anodes) 
  {
    num.regimes <- dim(regimes)[3]
    intervention.match <- array(dim = c(nrow(data), length(Anodes), 
                                        num.regimes))
    cum.intervention.match <- matrix(TRUE, nrow(data), num.regimes)
    for (Anode.index in seq_along(Anodes)) {
      cum.intervention.match <- cum.intervention.match & ((data[, 
                                                                Anodes[Anode.index]] == regimes[, Anode.index, ]) %in% 
                                                            c(TRUE, NA))
      intervention.match[, Anode.index, ] <- cum.intervention.match
    }
    return(intervention.match)
  }
CalcIPTW <-
  function (inputs, cum.g, msm.weights) 
  {
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", 
                    exact = TRUE))) {
      return(list(beta = NA, IC = matrix(NA, 1, 1)))
    }
    nodes <- inputs$all.nodes
    n <- nrow(inputs$data)
    num.regimes <- dim(inputs$regimes)[3]
    num.final.Ynodes <- length(inputs$final.Ynodes)
    Y.vec <- X.mat <- weight.vec <- NULL
    save.xy <- list()
    for (j in 1:num.final.Ynodes) {
      final.Ynode <- inputs$final.Ynodes[j]
      intervention.match <- InterventionMatch(inputs$intervention.match, 
                                              nodes$A, cur.node = final.Ynode)
      uncensored <- IsUncensored(inputs$uncensored, nodes$C, 
                                 final.Ynode)
      for (i in 1:num.regimes) {
        index <- uncensored & intervention.match[, i]
        col.index <- which.max(nodes$AC[nodes$AC < final.Ynode])
        Y <- inputs$data[index, final.Ynode]
        if (length(col.index > 0)) {
          g <- cum.g[index, col.index, i]
        }
        else {
          g <- 1
        }
        X <- inputs$combined.summary.measures[index, , i, 
                                              j]
        if (is.vector(X)) {
          dim(X) <- c(sum(index), ncol(inputs$combined.summary.measures))
        }
        weight <- msm.weights[index, i, j] * inputs$observation.weights[index]/g
        weight[msm.weights[index, i, j] == 0 | inputs$observation.weights[index] == 
                 0] <- 0
        save.xy[[length(save.xy) + 1]] <- list(X = X, Y = Y, 
                                               weight = weight, index = index)
        Y.vec <- c(Y.vec, Y)
        X.mat <- rbind(X.mat, X)
        weight.vec <- c(weight.vec, weight)
      }
    }
    colnames(X.mat) <- colnames(inputs$combined.summary.measures)
    if (nrow(X.mat) == 0) {
      warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
      num.beta <- ncol(inputs$combined.summary.measures)
      return(list(beta = rep(NA, num.beta), IC = matrix(nrow = n, 
                                                        ncol = num.beta)))
    }
    m.glm <- ltmle.glm(formula(inputs$working.msm), family = quasibinomial(), 
                       data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, 
                                                                                                  center = FALSE)))
    beta <- coef(m.glm)
    IC <- matrix(0, nrow = n, ncol = length(beta))
    m.beta <- array(dim = c(n, num.regimes, num.final.Ynodes))
    cnt <- 0
    for (j in 1:num.final.Ynodes) {
      final.Ynode <- inputs$final.Ynodes[j]
      for (i in 1:num.regimes) {
        newdata <- data.frame(inputs$combined.summary.measures[, 
                                                               , i, j])
        colnames(newdata) <- colnames(inputs$combined.summary.measures)
        SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, 
                                                        newdata = newdata, type = "response"), "prediction from a rank-deficient fit may be misleading")
        cnt <- cnt + 1
        XY.list <- save.xy[[cnt]]
        IC[XY.list$index, ] <- IC[XY.list$index, ] + XY.list$weight * 
          XY.list$X * (XY.list$Y - m.beta[XY.list$index, 
                                          i, j])
      }
    }
    C <- NormalizeIC(IC, inputs$combined.summary.measures, m.beta, 
                     msm.weights, observation.weights = inputs$observation.weights, 
                     g.ratio = NULL)
    normalized.IC <- t(safe.solve(C, t(IC)))
    household.IC <- HouseholdIC(normalized.IC, inputs$id)
    names(beta) <- inputs$beta.names
    return(list(beta = beta, IC = household.IC))
  }
CalcUncensoredMatrix <-
  function (data, Cnodes) 
  {
    uncensored <- matrix(nrow = nrow(data), ncol = length(Cnodes))
    cum.uncensored <- rep(TRUE, nrow(data))
    for (Cnode.index in seq_along(Cnodes)) {
      cum.uncensored <- cum.uncensored & (data[, Cnodes[Cnode.index]] %in% 
                                            c("uncensored", NA))
      uncensored[, Cnode.index] <- cum.uncensored
    }
    return(uncensored)
  }
CheckData <-
  function (data) 
  {
    if (inherits(data, "data.frame")) {
      return(droplevels(as.data.frame(data)))
    }
    stop("data must be a data frame (or inherit data.frame)")
  }
CheckForVarianceWarning <-
  function (inputs, g.ratio)
  {
    if (inputs$variance.method == "ic") {
      positivity <- mean(g.ratio < 1, na.rm = TRUE) > 0.01
      rare.events <- inputs$binaryOutcome && any(colMeans(inputs$data[,
                                                                      inputs$final.Ynodes, drop = FALSE], na.rm = TRUE) <
                                                   0.03)
      if (positivity || rare.events) {
        variance.available.warning <- VarianceAvailableWarning(inputs)
        warning.msg <- "Variance estimate is based on influence curve only, which may be significantly anticonservative because your data appears to contain"
        if (positivity)
          warning.msg <- paste(warning.msg, "positivity violations")
        if (positivity && rare.events)
          warning.msg <- paste(warning.msg, "and")
        if (rare.events)
          warning.msg <- paste(warning.msg, "rare events")
        if (is.null(variance.available.warning)) {
          warning.msg <- paste0(warning.msg, ". It is recommended to use variance.method='tmle' or variance.method='iptw' to obtain a more robust variance estimate (but run time may be significantly longer). See variance.method details in ?ltmle")
        }
        else {
          warning.msg <- paste0(warning.msg, ". ", variance.available.warning,
                                " but this will be addressed in a future release.")
        }
        warning(warning.msg)
      }
    }
    invisible(NULL)
  }
CheckInputs <-
  function (data, nodes, survivalOutcome, Qform, gform, gbounds,
            Yrange, deterministic.g.function, SL.library, SL.cvControl,
            regimes, working.msm, summary.measures, final.Ynodes, stratify,
            msm.weights, deterministic.Q.function, observation.weights,
            gcomp, variance.method, id){
    
    
    stopifnot(length(dim(regimes)) == 3)
    num.regimes <- dim(regimes)[3]
    if (!is.glm(GetLibrary(SL.library, "Q")) || !is.glm(GetLibrary(SL.library,
                                                                   "g"))) {
      if (!requireNamespace("SuperLearner"))
        stop("SuperLearner package is required if SL.library is not NULL or 'glm'")
    }
    if (!is.list(SL.cvControl)) {
      stop("SL.cvControl must be a list")
    }
    ## if (!(all(names(SL.cvControl) %in% c("V", "stratifyCV", "shuffle")))) {
    ## stop("The valid names for SL.cvControl are V, stratifyCV, shuffle. validRows is not currently supported. See ?SuperLearner.CV.control")
    ## }
    if (is.unsorted(nodes$A, strictly = TRUE))
      stop("Anodes must be in increasing order")
    if (is.unsorted(nodes$C, strictly = TRUE))
      stop("Cnodes must be in increasing order")
    if (is.unsorted(nodes$L, strictly = TRUE)){
      lnodes = names(data)[nodes$L]
      lpos = nodes$L
      print(rbind(lnodes,lpos))
      stop("Lnodes must be in increasing order. See wrong order above.\n")
    }
    if (is.unsorted(nodes$Y, strictly = TRUE))
      stop("Ynodes must be in increasing order")
    if (is.unsorted(final.Ynodes, strictly = TRUE))
      stop("final.Ynodes must be in increasing order")
    if (length(nodes$L) > 0) {
      if (max(nodes$L) > max(nodes$Y))
        stop("Lnodes are not allowed after the final Y node")
    }
    all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
    if (length(all.nodes) > length(unique(all.nodes)))
      stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes",
           "\nThe following variables are listed in more than one Xnodes list:\n",
           paste0(names(data)[all.nodes[duplicated(all.nodes)]],collapse=", "))
    if (is.null(nodes$Y))
      stop("Ynodes cannot be null")
    if (length(nodes$AC) > 0 && !all(sseq(min(nodes$AC), ncol(data)) %in%
                                     all.nodes)) {
      tmp=sseq(min(nodes$AC), ncol(data))
      stop(paste0("The following variables are not listed as A-, C-, L-, or Y-nodes:\n",
                  paste(names(data)[tmp[!(tmp%in%all.nodes)]],collapse=", "),
                  "\n",
                  "All nodes after the first A/C node must be in A-, C-, L-, or Ynodes"))
      ## stop("All nodes after the first A/C node must be in A-, C-, L-, or Ynodes")
    }
    for (reserved.name in c("observation.weights", "Q.kplus1",
                            "Qstar")) {
      if (reserved.name %in% names(data))
        stop(paste(reserved.name, "is reserved and may not be used as a column name of data"))
    }
    if (length(variance.method) != 1 || !(variance.method %in%
                                          c("ic", "tmle", "iptw"))) {
      stop("variance.method must be one of 'ic', 'tmle', 'iptw'")
    }
    if (length(gform) > 0) {
      if (is.character(gform)) {
        if (length(gform) != length(nodes$AC)){
          stop(paste0("length(gform) != length(c(Anodes, Cnodes))",
                      "\n   gform: ",
                      paste0(gform,collapse = "\n          "),
                      "\nAC nodes: ",
                      paste0(names(data)[nodes$AC],collapse = ", ")))
        }
        for (i in 1:length(gform)) {
          if (LhsVars(gform[i]) != names(data)[nodes$AC[i]]) {
            stop("The LHS of gform[", i, "] should be the name of the ",
                 i, "th A or C node.")
          }
          parents <- if (nodes$AC[i] > 1) {
            names(data)[1:(nodes$AC[i] - 1)]
          }
          else {
            NULL
          }
          if (!all(RhsVars(gform[i]) %in% parents)) {
            stop("Some nodes in gform[", i, "] are not parents of ",
                 LhsVars(gform[i]))
          }
          if (any(RhsVars(gform[i]) %in% names(data)[nodes$C]))
            stop("Cnodes should not be used as RHS variables in gform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
        }
      }
      else {
        if (!is.numeric(gform))
          stop("gform should be a character vector or numeric")
        if (nrow(gform) != nrow(data))
          stop("if gform is numeric, it should have the same number of rows as data")
        if (ncol(gform) != length(nodes$AC))
          stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
        if (length(dim(gform)) != 3 || dim(gform)[3] != num.regimes)
          stop("if gform is numeric, dim[3] should be num.regimes (gform can also be a matrix if variance.method == 'ic')")
        if (!is.null(deterministic.g.function))
          stop("if gform is numeric, deterministic.g.function must be NULL")
        if (max(gform, na.rm = T) > 1 || min(gform, na.rm = T) <
            0)
          stop("if gform is numeric, all values should be probabilities")
        if (!is.null(deterministic.Q.function) && !isTRUE(attr(data,
                                                               "called.from.estimate.variance", exact = TRUE)))
          warning("If gform is numeric and deterministic.Q.function is not NULL, deterministic.Q.function will only affect g based on the observed values of the Anodes, not the counterfactual values. If your deterministic.Q.function does depends on the values of the Anodes, it is recommended to not use numeric gform.")
      }
    }
    if (length(Qform) > 0) {
      if (!is.character(Qform))
        stop("Qform should be a character vector")
      if (length(Qform) != length(nodes$LY))
        stop(paste0("length of Qform is not equal to number of L/Y nodes\n",
                    "  Qforms: ",paste0(names(Qform),collapse = ", "),"\n",
                    "LY nodes: ",paste0(names(data)[nodes$LY],collapse = ", ")))
      for (i in 1:length(Qform)) {
        if (length(names(Qform[i])) == 0)
          stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
        if (names(Qform[i]) != names(data)[nodes$LY[i]])
          stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
        if (Qform[i] != "IDENTITY") {
          if (LhsVars(Qform[i]) != "Q.kplus1")
            stop("LHS of each Qform should be Q.kplus1")
          parents <- names(data)[1:(nodes$LY[i] - 1)]
          if (!all(RhsVars(Qform[i]) %in% parents)) {
            stop("Some nodes in Qform[", i, "] are not parents of ",
                 names(Qform[i]))
          }
          if (any(RhsVars(Qform[i]) %in% names(data)[nodes$C]))
            stop("Cnodes should not be used as RHS variables in Qform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
        }
      }
    }
    if (length(gbounds) != 2)
      stop("gbounds should have length 2")
    if (!(is.null(deterministic.g.function) || is.function(deterministic.g.function)))
      stop("deterministic.g.function should be a function or NULL")
    if (!all(unlist(data[, nodes$A]) %in% c(0, 1, NA)))
      stop("in data, all Anodes should be binary")
    if (!all(sapply(data[, c(nodes$A, nodes$Y), drop = F], is.numeric)))
      stop("in data, all Anodes and Ynodes should be numeric (not, for instance, logical)")
    if (any(sapply(data, is.infinite)))
      stop("infinite values are not supported in data")
    all.Y <- unlist(data[, nodes$Y])
    binaryOutcome <- all(all.Y %in% c(0, 1, NA))
    if (binaryOutcome) {
      if (is.null(survivalOutcome)) {
        if (length(nodes$Y) == 1) {
          survivalOutcome <- FALSE
        }
        else {
          stop("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified if there are multiple Ynodes.")
        }
      }
      if (!is.null(Yrange) && !is.equal(Yrange, c(0L, 1L))) {
        stop("All Ynodes are 0, 1, or NA, but Yrange is something other than NULL or c(0, 1)")
      }
    }
    else {
      if (is.null(survivalOutcome))
        survivalOutcome <- FALSE
      if (survivalOutcome)
        stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
    }
    uncensored.array <- CalcUncensoredMatrix(data, nodes$C)
    for (Ynode in nodes$Y) {
      uncensored <- IsUncensored(uncensored.array, nodes$C,
                                 cur.node = Ynode)
      deterministic <- IsDeterministic(data, cur.node = Ynode,
                                       deterministic.Q.function = NULL, nodes, called.from.estimate.g = FALSE,
                                       survivalOutcome)$is.deterministic
      if (anyNA(data[deterministic, Ynode]) || !all(data[deterministic,
                                                         Ynode] == 1))
        stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")
      if (anyNA(data[uncensored, Ynode]))
        stop("Ynodes may not be NA except after censoring")
    }
    if (!is.equal(dim(regimes)[1:2], c(nrow(data), length(nodes$A))))
      stop("Problem with abar or regimes:\n   In ltmleMSM, regimes should have dimensions n x num.Anodes x num.regimes\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
    stopifnot(num.regimes == nrow(summary.measures))
    if (!all(regimes %in% c(0, 1, NA)))
      stop("all regimes should be binary")
    for (Anode.index in seq_along(nodes$A)) {
      first.LYnode <- min(nodes$LY)
      cur.node <- nodes$A[Anode.index]
      uncensored <- IsUncensored(uncensored.array, Cnodes = nodes$C,
                                 cur.node = cur.node)
      deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function,
                                       nodes, called.from.estimate.g = TRUE, survivalOutcome)$is.deterministic
      if (anyNA(regimes[uncensored & !deterministic, Anode.index,
      ])) {
        stop("NA in regimes/abar not allowed (except after censoring/death)")
      }
      if (cur.node < first.LYnode && anyNA(regimes[!deterministic,
                                                   Anode.index, ])) {
        warning("NA in regimes/abar before the first L/Y node will probably cause an error")
      }
    }
    num.final.Ynodes <- length(final.Ynodes)
    if ((length(dim(summary.measures)) != 3) || !is.equal(dim(summary.measures)[c(1,
                                                                                  3)], c(num.regimes, num.final.Ynodes)))
      stop("summary.measures should be an array with dimensions num.regimes x num.summary.measures x num.final.Ynodes")
    if (class(working.msm) != "character")
      stop("class(working.msm) must be 'character'")
    if (LhsVars(working.msm) != "Y")
      stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
    if (!is.vector(observation.weights) || length(observation.weights) !=
        nrow(data) || anyNA(observation.weights) || any(observation.weights <
                                                        0) || max(observation.weights) == 0)
      stop("observation.weights must be NULL or a vector of length nrow(data) with no NAs, no negative values, and at least one positive value")
    if (!(is.null(msm.weights) || is.equal(msm.weights, "empirical") ||
          is.equal(dim(msm.weights), c(nrow(data), num.regimes,
                                       num.final.Ynodes)) || is.equal(dim(msm.weights),
                                                                      c(num.regimes, num.final.Ynodes)))) {
      stop("msm.weights must be NULL, 'empirical', or an array with dim(msm.weights) = c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes)")
    }
    if (!(is.null(id) || ((is.factor(id) || is.vector(id)) &&
                          length(id) == nrow(data))))
      stop("id must be a vector with length nrow(data) or be NULL")
    return(list(survivalOutcome = survivalOutcome, binaryOutcome = binaryOutcome,
                uncensored = uncensored.array))
  }
CheckVarianceEstimateRatio <-
  function (summary.obj) 
  {
    if (anyNA(summary.obj$variance.estimate.ratio)) {
      warning("Unable to compute standard errors.")
      return(NULL)
    }
    if (any(summary.obj$variance.estimate.ratio > 100)) {
      warning.msg <- paste0("max(TMLE based variance estimate / IC based variance estimate) = ", 
                            floor(max(summary.obj$variance.estimate.ratio)), 
                            ".\nWhen this ratio is greater than 100, both variance estimates are less likely to be accurate.")
      warning(warning.msg)
    }
  }
CleanData <- function(data,
                      nodes,
                      deterministic.Q.function,
                      survivalOutcome,
                      showMessage = TRUE){
  is.nan.df <- function(x) {
    y <- if (length(x)) {
      do.call("cbind", lapply(x, "is.nan"))
    }
    else {
      matrix(FALSE, length(row.names(x)), 0)
    }
  }
  is.na.strict <- function(x) is.na(x) & !is.nan.df(x)
  changed <- FALSE
  ua <- rep(TRUE, nrow(data))
  if (ncol(data) == 1)
    return(data)
  deterministic.Q.function.depends.on.called.from.estimate.g <- !is.null(deterministic.Q.function) &&
    length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) >
    0
  for (i in 1:(ncol(data) - 1)) {
    if (anyNA(data[ua, i])){
      stop(paste0("Missing values in variable ",names(data)[i],".\n",
                  "NA values are not permitted in data except after censoring or a survival event"))
    }
    ## stop("NA values are not permitted in data except after censoring or a survival event")
    if (i %in% c(nodes$L, nodes$Y, nodes$AC)) {
      is.deterministic <- ua & IsDeterministic(data, cur.node = i +
                                                 1, deterministic.Q.function = deterministic.Q.function,
                                               nodes = nodes, called.from.estimate.g = TRUE,
                                               survivalOutcome = survivalOutcome)$is.deterministic
      if (deterministic.Q.function.depends.on.called.from.estimate.g) {
        is.deterministic.Q <- ua & IsDeterministic(data,
                                                   cur.node = i + 1, deterministic.Q.function = deterministic.Q.function,
                                                   nodes = nodes, called.from.estimate.g = FALSE,
                                                   survivalOutcome = survivalOutcome)$is.deterministic
        if (any(is.deterministic[ua] & !is.deterministic.Q[ua]))
          stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)")
      }
      ua[ua] <- !is.deterministic[ua]
      if (anyNA(ua))
        stop("internal ltmle error - ua should not be NA in CleanData")
      if (!all(is.na.strict(data[is.deterministic, setdiff((i +
                                                            1):ncol(data), nodes$Y), drop = FALSE]))) {
        data[is.deterministic, setdiff((i + 1):ncol(data),
                                       nodes$Y)] <- NA
        changed <- TRUE
      }
      if (i %in% nodes$C) {
        censored <- data[, i] == "censored" & ua
        if (!all(is.na.strict(data[censored, (i + 1):ncol(data),
                                   drop = FALSE]))) {
          data[censored, (i + 1):ncol(data)] <- NA
          changed <- TRUE
        }
        ua[ua] <- !censored[ua]
        if (anyNA(ua))
          stop("internal ltmle error - ua should not be NA in CleanData")
      }
      if (changed && showMessage) {
        message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging ltmle.")
      }
    }
  }
  return(data)
}
ConvertCensoringNodes <-
  function (data, Cnodes, has.deterministic.functions = FALSE) 
  {
    error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored'\n See ?BinaryToCensoring \n (binary is also accepted, where 0=censored, 1=uncensored, but this is not recommended)"
    for (i in Cnodes) {
      col <- data[, i]
      if (is.numeric(col)) {
        if (!all(col %in% c(0, 1, NA))) 
          stop(error.msg)
        data[, i] <- BinaryToCensoring(is.uncensored = col)
        if (has.deterministic.functions) 
          warning("Censoring nodes have been converted from binaries to factors - see ?BinaryToCensoring.\n Note that if you are writing your own deterministic.g.function or deterministic.Q.function that censoring nodes are converted to factors\n before these functions are called.")
      }
      else if (is.factor(col)) {
        if (!all(levels(col) %in% c("censored", "uncensored"))) {
          stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
        }
      }
      else {
        stop(error.msg)
      }
    }
    return(data)
  }
ConvertCensoringNodeToBinary <-
  function (x) 
  {
    stopifnot(is.factor(x) && all(levels(x) %in% c("censored", 
                                                   "uncensored")))
    b <- rep(NA_integer_, length(x))
    b[x == "censored"] <- 0L
    b[x == "uncensored"] <- 1L
    return(b)
  }
ConvertToMainTerms <-
  function (data, msm, summary.measures, nodes) 
  {
    baseline.column.names <- names(data)[nodes$baseline]
    summary.column.names <- colnames(summary.measures)
    rhs.vars <- RhsVars(msm)
    if (length(intersect(baseline.column.names, summary.column.names)) > 
        0) 
      stop("Baseline covariate columns of data and columns of summary.measures may not have the same name")
    if (!all(rhs.vars %in% c(baseline.column.names, summary.column.names))) 
      stop("All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
    baseline.column.names <- intersect(baseline.column.names, 
                                       rhs.vars)
    baseline.data <- data[, baseline.column.names, drop = FALSE]
    num.regimes <- dim(summary.measures)[1]
    num.summary.measures <- dim(summary.measures)[2]
    num.final.Ynodes <- dim(summary.measures)[3]
    n <- nrow(data)
    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimes) {
        combined.summary.measures <- model.matrix(as.formula(msm), 
                                                  data.frame(Y = 1, baseline.data, matrix(summary.measures[i, 
                                                                                                           , j], nrow = n, ncol = num.summary.measures, 
                                                                                          byrow = TRUE, dimnames = list(NULL, colnames(summary.measures)))))
        if (i == 1 && j == 1) {
          main.terms.summary.measures <- array(dim = c(n, 
                                                       ncol(combined.summary.measures), num.regimes, 
                                                       num.final.Ynodes))
          beta.names <- colnames(combined.summary.measures)
        }
        main.terms.summary.measures[, , i, j] <- combined.summary.measures
      }
    }
    colnames(main.terms.summary.measures) <- paste("S", 1:ncol(main.terms.summary.measures), 
                                                   sep = "")
    main.terms.msm <- paste("Y ~ -1 +", paste(colnames(main.terms.summary.measures), 
                                              collapse = " + "))
    return(list(msm = main.terms.msm, summary.measures = main.terms.summary.measures, 
                beta.names = beta.names, baseline.column.names = baseline.column.names))
  }
CreateInputs <- function (data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome,
                          Qform, gform, gbounds, Yrange, deterministic.g.function,
                          SL.library, SL.cvControl, regimes, working.msm, summary.measures,
                          final.Ynodes, stratify, msm.weights, estimate.time, gcomp,
                          iptw.only, deterministic.Q.function, variance.method, observation.weights,
                          id,verbose){
  
  
  
  if (is.list(regimes)) {
    if (!all(sapply(regimes, is.function)))
      stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- simplify2array(lapply(regimes, function(rule) drop3(RegimesFromAbar(data,
                                                                                   rule = rule))), higher = TRUE)
  }
  if (!(is.null(regimes) || length(dim(regimes)) == 3)) {
    stop("regimes must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")
  }
  if (is.null(regimes) || dim(regimes)[3] == 0) {
    if (length(Anodes) != 0) {
      stop("regimes must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")
    }
    regimes <- array(numeric(0), dim = c(nrow(data), 0, 1))
  }
  num.regimes <- dim(regimes)[3]
  if (is.logical(regimes)) {
    regimes <- regimes * 1
    message("abar or regimes was passed as logical and was converted to numeric")
  }
  #browser()
  all.nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, all.nodes, check.Qform = TRUE, Qform = Qform)$Qform
  
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions = !is.null(deterministic.g.function) &&
                                  is.null(deterministic.Q.function))
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(all.nodes$Y)
  }
  else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
  if (identical(SL.library, "default"))
    SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")
  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow = num.regimes, ncol = 0)
  }
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m = 1,
                                     n = num.final.Ynodes), dim = c(nrow(summary.measures),
                                                                    ncol(summary.measures), num.final.Ynodes), dimnames = list(rownames(summary.measures),
                                                                                                                               colnames(summary.measures), NULL))
  }
  if (is.null(observation.weights))
    observation.weights <- rep(1, nrow(data))
  if (is.matrix(gform)) {
    if (num.regimes > 1 && variance.method != "ic")
      stop("If there is more than one regime (using ltmle with list abar or ltmleMSM) and numeric gform and variance.method != 'ic', then gform must be an array, not a matrix.")
    gform <- array(gform, dim = c(nrow(gform), ncol(gform),
                                  num.regimes))
  }
  check.results <- CheckInputs(data, all.nodes, survivalOutcome,
                               Qform, gform, gbounds, Yrange, deterministic.g.function,
                               SL.library, SL.cvControl, regimes, working.msm, summary.measures,
                               final.Ynodes, stratify, msm.weights, deterministic.Q.function,
                               observation.weights, gcomp, variance.method, id)
  survivalOutcome <- check.results$survivalOutcome
  if (!isTRUE(attr(data, "called.from.estimate.variance", exact = TRUE)) &&
      !isTRUE(attr(data, "skip.clean.data", exact = TRUE))) {
    data <- CleanData(data, all.nodes, deterministic.Q.function,
                      survivalOutcome)
  }
  transform.list <- TransformOutcomes(data, all.nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  binaryOutcome <- check.results$binaryOutcome
  if (length(Qform) == 0)
    Qform <- GetDefaultForm(data, all.nodes, is.Qform = TRUE,
                            stratify, survivalOutcome, showMessage = TRUE)
  if (length(gform) == 0)
    gform <- GetDefaultForm(data, all.nodes, is.Qform = FALSE,
                            stratify, survivalOutcome, showMessage = TRUE)
  main.terms <- ConvertToMainTerms(data, working.msm, summary.measures,
                                   all.nodes)
  intervention.match <- CalcInterventionMatchArray(data, regimes,
                                                   all.nodes$A)
  inputs <- list(data = data, all.nodes = all.nodes, survivalOutcome = survivalOutcome,
                 Qform = Qform, gform = gform, gbounds = gbounds, Yrange = Yrange,
                 deterministic.g.function = deterministic.g.function,
                 SL.library.Q = SL.library.Q, SL.library.g = SL.library.g,
                 SL.cvControl = SL.cvControl, regimes = regimes, working.msm = main.terms$msm,
                 combined.summary.measures = main.terms$summary.measures,
                 final.Ynodes = final.Ynodes, stratify = stratify, msm.weights = msm.weights,
                 estimate.time = estimate.time, gcomp = gcomp, iptw.only = iptw.only,
                 deterministic.Q.function = deterministic.Q.function,
                 binaryOutcome = binaryOutcome, transformOutcome = transformOutcome,
                 variance.method = variance.method, observation.weights = observation.weights,
                 baseline.column.names = main.terms$baseline.column.names,
                 beta.names = main.terms$beta.names, uncensored = check.results$uncensored,
                 intervention.match = intervention.match, id = id,verbose=verbose)
  class(inputs) <- "ltmleInputs"
  if (length(all.nodes$AC) == 0)
    inputs$variance.method <- "ic"
  if (inputs$variance.method != "ic" && !is.null(VarianceAvailableWarning(inputs)))
    inputs$variance.method <- "ic"
  
  
  return(inputs)
}
CreateLYNodes <-
  function (data, nodes, check.Qform, Qform) {
    
    
    LYnodes <- sort(c(nodes$L, nodes$Y))
    SuppressGivenWarnings(nodes.to.remove <- LYnodes[LYnodes < 
                                                       min(nodes$AC)], "no non-missing arguments to min; returning Inf")
    if(length(LYnodes) > 1){
      for(i in 1:(length(LYnodes) - 1)){
        cur.node <- LYnodes[i]
        next.node <- LYnodes[i + 1]
        if(!any(cur.node:next.node %in% nodes$AC)){
          nodes.to.remove <- c(nodes.to.remove, next.node)
        }
      }
    }
    new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
    if (check.Qform) {
      removed.Qform.index <- NULL
      for(i in nodes.to.remove){
        index <- which(names(Qform) == names(data)[i])
        if(length(index) > 0){
          removed.Qform.index <- c(removed.Qform.index, 
                                   index)
        }
      }
      if (!is.null(removed.Qform.index)) {
        message("L/Y nodes (after removing blocks)  : ", 
                paste(names(data)[new.LYnodes], collapse = " "), 
                "\n")
        message("Qform names                        : ", 
                paste(names(Qform), collapse = " "), "\n")
        message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block\nof L/Y nodes. They are being dropped from Qform:\n"), 
                paste(names(Qform)[removed.Qform.index], "\n", 
                      collapse = " "))
        Qform <- Qform[-removed.Qform.index]
      }
      return(list(LYnodes = new.LYnodes, Qform = Qform))
    }
    return(new.LYnodes)
  }
CreateNodes <-
  function (data, Anodes, Cnodes, Lnodes, Ynodes) 
  {
    Anodes <- NodeToIndex(data, Anodes)
    Cnodes <- NodeToIndex(data, Cnodes)
    Lnodes <- NodeToIndex(data, Lnodes)
    Ynodes <- NodeToIndex(data, Ynodes)
    nodes <- SuppressGivenWarnings(list(A = Anodes, C = Cnodes, 
                                        L = Lnodes, Y = Ynodes, AC = sort(c(Anodes, Cnodes))), 
                                   "is.na() applied to non-(list or vector) of type 'NULL'")
    nodes$baseline <- sseq(1, min(c(nodes$A, nodes$L, nodes$C, 
                                    nodes$Y)) - 1)
    
    #browser()
    nodes$LY <- CreateLYNodes(data, nodes, check.Qform = FALSE)
    return(nodes)
  }
deterministic.g.function_template <-
  function (data, current.node, nodes) 
  {
    is.deterministic <- stop("replace me!")
    prob1 <- stop("replace me!")
    return(list(is.deterministic = is.deterministic, prob1 = prob1))
  }
deterministic.Q.function_template <-
  function (data, current.node, nodes, called.from.estimate.g) 
  {
    is.deterministic <- stop("replace me!")
    Q.value <- stop("replace me!")
    return(list(is.deterministic = is.deterministic, Q.value = Q.value))
  }
drop3 <-
  function (x) 
  {
    return(dropn(x, 3))
  }
dropn <-
  function (x, n) 
  {
    stopifnot(length(dim(x)) == n)
    stopifnot(dim(x)[n] == 1)
    dn <- dimnames(x)
    dim(x) <- dim(x)[1:(n - 1)]
    dimnames(x) <- dn[1:(n - 1)]
    return(x)
  }
# called from EstimateG
Estimate <- function(inputs,
                     form,
                     subs,
                     family,
                     type,
                     nodes,
                     Qstar.kplus1,
                     cur.node,
                     calc.meanL,
                     called.from.estimate.g,
                     regimes.meanL,
                     regimes.with.positive.weight)
{
  FitAndPredict <- function() {
    if (length(Y.subset) < 2)
      stop("Estimation failed because there are fewer than 2 observations to fit")
    Y.subset.range <- range(Y.subset)
    if (anyNA(Y.subset.range))
      stop("Internal error - NA in Y during Estimate")
    if (Y.subset.range[1] < -0.0001 || Y.subset.range[2] >
        1.0001)
      stop("Internal error - Y negative or greater than 1 in Estimate")
    if (Y.subset.range[2] - Y.subset.range[1] < 0.0001) {
      Y.value <- Y.subset.range[2]
      m <- list("no estimation occured because all Y values are the same",
                Y.value = Y.value)
      predicted.values <- ValuesByType(rep(Y.value, nrow(newdata)))
      class(m) <- "no.Y.variation"
    }
    else {
      if (use.glm) {
        if (inputs$verbose){ message("Estimate: calling ltmle.glm.fit ...")}
        SuppressGivenWarnings({
          m <- ltmle.glm.fit(y = Y.subset, x = X.subset,
                             family = family, weights = observation.weights.subset,
                             offset = offst, intercept = intercept)
          ## if (inputs$verbose){ message("Fitted coefs: ")
          ## browser()
          ## print(paste("Sample size: ",m$n))
          ## print(coef(m))
          ## }
          m$terms <- tf
          predicted.values <- predict(m,newdata = newdata,type = type)
        }, GetWarningsToSuppress())
      }
      else {
        newX.list <- GetNewX(newdata)
        SetSeedIfRegressionTesting()
        if (SL.library[[1]]=="glmnet"){
          if (inputs$verbose){ message("Estimate: calling ltmle.glmnet ...")}
          try.result <- try({
            m <- ltmle.glmnet(Y=Y.subset,
                              X=X.subset,
                              newX=newX.list$newX,
                              family=family,
                              obsWeights=observation.weights.subset,
                              id =id.subset,
                              alpha=inputs$SL.cvControl$alpha,
                              selector=inputs$SL.cvControl$selector,
            )
          })
          predicted.values <- ProcessSLPrediction(pred=m$predicted.values,
                                                  new.subs=newX.list$new.subs,
                                                  try.result=try.result)
          m <- m$fit
        }else{
          try.result <- try({
            SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y = Y.subset,
                                                                  X = X.subset, SL.library = SL.library, cvControl = inputs$SL.cvControl,
                                                                  verbose = FALSE, family = family, newX = newX.list$newX,
                                                                  obsWeights = observation.weights.subset,
                                                                  id = id.subset, env = environment(SuperLearner::SuperLearner)),
                                  c("non-integer #successes in a binomial glm!",
                                    "prediction from a rank-deficient fit may be misleading"))
          })
          predicted.values <- ProcessSLPrediction(m$SL.predict,
                                                  newX.list$new.subs, try.result)
        }
      }
    }
    return(list(m = m, predicted.values = predicted.values))
  }
  GetSLStopMsg <- function(Y) {
    ifelse(all(Y %in% c(0, 1, NA)), "", "\n Note that some SuperLeaner libraries crash when called with continuous dependent variables, as in the case of initial Q regressions with continuous Y or subsequent Q regressions even if Y is binary.")
  }
  ProcessSLPrediction <- function(pred, new.subs, try.result) {
    if (inherits(try.result, "try-error")) {
      stop(paste("\n\nError occured during call to SuperLearner:\n",
                 form, GetSLStopMsg(Y.subset), "\n The error reported is:\n",
                 try.result))
    }
    if (all(is.na(pred))) {
      stop(paste("\n\n Unexpected error: SuperLearner returned all NAs during regression:\n",
                 form, GetSLStopMsg(Y.subset)))
    }
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- pred
    if (max(predicted.values, na.rm = T) > 1 || min(predicted.values,
                                                    na.rm = T) < 0) {
      msg <- paste("SuperLearner returned predicted.values > 1 or < 0: [min, max] = [",
                   min(predicted.values, na.rm = T), ",", max(predicted.values,
                                                              na.rm = T), "]. Bounding to [0,1]")
      warning(msg)
      predicted.values <- Bound(predicted.values, bounds = c(0,
                                                             1))
    }
    return(ValuesByType(predicted.values))
  }
  PredictOnly <- function(newdata1) {
    if (class(m)[1] == "no.Y.variation")
      return(rep(m$Y.value, nrow(newdata1)))
    if (use.glm) {
      SuppressGivenWarnings(pred <- predict(m, newdata1,
                                            type), "prediction from a rank-deficient fit may be misleading")
    }
    else {
      if  (SL.library[[1]]=="glmnet"){
        newX.list <- GetNewX(newdata1)
        pred <- ProcessSLPrediction(pred=predict(m,newX=newX.list$newX),
                                    new.subs=newX.list$new.subs,
                                    try.result = NULL)
        
      }else{
        newX.list <- GetNewX(newdata1)
        pred <- ProcessSLPrediction(predict(m,
                                            newX.list$newX,
                                            X.subset,
                                            Y.subset,
                                            onlySL = TRUE)$pred,
                                    newX.list$new.subs,
                                    try.result = NULL)
      }
    }
    return(pred)
  }
  ValuesByType <- function(x) {
    if (type == "link") {
      stopifnot(family$family %in% c("binomial", "quasibinomial"))
      qlogis(Bound(x, bounds = c(0.0001, 0.9999)))
    }
    else {
      x
    }
  }
  GetNewX <- function(newdata1) {
    new.mod.frame <- model.frame(f, data = newdata1, drop.unused.levels = TRUE,
                                 na.action = na.pass)
    newX.temp <- model.matrix(terms(f), new.mod.frame)
    Xnames=colnames(newX.temp)
    if (!use.glm) {
      colnames(newX.temp) <- paste0("Xx.", 1:ncol(newX.temp))
    }
    new.subs <- !matrixStats::rowAnyMissings(newX.temp)
    newX <- as.data.frame(newX.temp[new.subs, , drop = FALSE])
    if (ncol(X) == 1) {
      X.subset <<- cbind(X.subset, ltmle.added.constant = 1)
      newX <- cbind(newX, ltmle.added.constant = 1)
    }
    attr(newX,"Xnames") <- Xnames
    return(list(newX = newX, new.subs = new.subs))
  }
  PredictProbAMeanL <- function() {
    probAis1.meanL <- matrix(NaN, nrow(inputs$data), length(nodes$LY) -
                               1)
    if (ncol(probAis1.meanL) == 0)
      return(probAis1.meanL)
    all.LY.nodes <- sort(union(nodes$L, nodes$Y))
    LYindex <- length(nodes$LY)
    for (i in length(all.LY.nodes):1) {
      regression.node <- all.LY.nodes[i]
      L <- data[single.subs, regression.node]
      if (is.numeric(L) && !IsBinary(L)) {
        meanL <- mean(L, na.rm = TRUE)
      }
      else {
        meanL <- Mode(L, na.rm = TRUE)
      }
      newdata.meanL[, regression.node] <- meanL
      if (regression.node %in% nodes$LY[1:length(nodes$LY) -
                                        1]) {
        LYindex <- LYindex - 1
        probAis1.meanL[, LYindex] <- PredictOnly(newdata = newdata.meanL)
      }
    }
    if (anyNA(probAis1.meanL[, 1]))
      stop("NA in probAis1.meanL[, 1]")
    return(probAis1.meanL)
  }
  #
  # Body of function Estimate starts here
  #
  stopifnot(type %in% c("link", "response"))
  num.regimes <- dim(inputs$regimes)[3]
  if (form == "IDENTITY") {
    stopifnot(is.vector(Qstar.kplus1) == 1)
    predicted.values <- ValuesByType(matrix(Qstar.kplus1,
                                            nrow = nrow(inputs$data), ncol = num.regimes))
    fit <- as.list(rep("no fit because form == IDENTITY",
                       num.regimes))
    deterministic.list.olddata <- IsDeterministic(inputs$data,
                                                  cur.node, inputs$deterministic.Q.function, nodes,
                                                  called.from.estimate.g, inputs$survivalOutcome)
    is.deterministic <- matrix(deterministic.list.olddata$is.deterministic,
                               nrow = nrow(inputs$data), ncol = num.regimes)
    deterministic.Q <- matrix(NA, nrow(inputs$data), num.regimes)
    deterministic.Q[is.deterministic, ] <- deterministic.list.olddata$Q
    return(list(predicted.values = predicted.values, fit = fit,
                is.deterministic = is.deterministic, deterministic.Q = deterministic.Q,
                prob.A.is.1.meanL = NULL))
  }
  data <- inputs$data
  if (cur.node %in% nodes$C) {
    data[, cur.node] <- ConvertCensoringNodeToBinary(data[, cur.node])
  }
  f <- as.formula(form)
  SL.library <- if (called.from.estimate.g)
    inputs$SL.library.g
  else inputs$SL.library.Q
  use.glm <- (is.glm(SL.library) || length(RhsVars(f)) == 0)
  first.regime <- min(regimes.with.positive.weight)
  if (is.null(Qstar.kplus1)) {
    data.with.Qstar <- data
  }
  else {
    if (is.matrix(Qstar.kplus1)) {
      data.with.Qstar <- cbind(data, Q.kplus1 = Qstar.kplus1[,
                                                             first.regime])
    }
    else {
      data.with.Qstar <- cbind(data, Q.kplus1 = Qstar.kplus1)
    }
  }
  if (inputs$verbose){ message("Estimate: framing formula ",form," into Y and X...")}
  mod.frame <- model.frame(f, data = data.with.Qstar, drop.unused.levels = TRUE,
                           na.action = na.pass)
  Y <- mod.frame[[1]]
  tf <- terms(f)
  X <- model.matrix(tf, mod.frame)
  offst <- model.offset(mod.frame)
  intercept <- attributes(tf)$intercept
  if (!use.glm) {
    if (is.equal(family, quasibinomial()))
      family <- binomial()
    if (!is.null(offst))
      stop("offset in formula not supported with SuperLearner")
    colnames(X) <- paste0("Xx.", 1:ncol(X))
    X <- as.data.frame(X)
  }
  fit <- vector("list", num.regimes)
  predicted.values <- deterministic.Q <- matrix(NA, nrow(data),
                                                num.regimes)
  is.deterministic <- matrix(FALSE, nrow(data), num.regimes)
  fit.and.predict <- NULL
  multiple.subs <- is.matrix(subs)
  multiple.Qstar <- is.matrix(Qstar.kplus1)
  if (calc.meanL) {
    prob.A.is.1.meanL <- array(NaN, dim = c(nrow(inputs$data),
                                            num.regimes, length(nodes$LY) - 1))
    Anode.index <- which(nodes$A < cur.node)
  }
  else {
    prob.A.is.1.meanL <- NULL
  }
  for (regime.index in regimes.with.positive.weight) {
    newdata <- SetA(data = data.with.Qstar, regimes = inputs$regimes,
                    Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)
    if (calc.meanL) {
      if (!is.null(regimes.meanL)) {
        newdata.meanL <- SetA(data = data.with.Qstar,
                              regimes = regimes.meanL, Anodes = nodes$A,
                              regime.index = regime.index, cur.node = cur.node)
      }
      else {
        newdata.meanL <- newdata
      }
    }
    deterministic.list.newdata <- IsDeterministic(newdata,
                                                  cur.node, inputs$deterministic.Q.function, nodes,
                                                  called.from.estimate.g, inputs$survivalOutcome)
    if (called.from.estimate.g && !is.null(inputs$deterministic.g.function)) {
      newdata.with.current <- newdata
      stopifnot(cur.node %in% nodes$AC)
      if (cur.node %in% nodes$A) {
        newdata.with.current[, cur.node] <- inputs$regimes[,
                                                           which(nodes$A == cur.node), regime.index]
      }
      else {
        newdata.with.current <- newdata
      }
      deterministic.g.list.newdata <- IsDeterministicG(newdata.with.current,
                                                       cur.node, inputs$deterministic.g.function, nodes,
                                                       using.newdata = T)
    }
    else {
      deterministic.g.list.newdata <- list(is.deterministic = rep(FALSE,
                                                                  nrow(data)), prob1 = NULL)
    }
    if (regime.index > first.regime && multiple.Qstar) {
      Y <- Qstar.kplus1[, regime.index]
    }
    if (regime.index == first.regime || multiple.subs) {
      single.subs <- if (multiple.subs)
        subs[, regime.index]
      else subs
      X.subset <- X[single.subs, , drop = FALSE]
      id.subset <- inputs$id[single.subs]
      ## if (any(is.na(single.subs))) browser()
      if (any(single.subs))
        X.subset[, colAlls(X.subset == 0)] <- 1
      observation.weights.subset <- inputs$observation.weights[single.subs]
      offst.subset <- offst[single.subs]
    }
    if (regime.index == first.regime || multiple.subs ||
        multiple.Qstar) {
      Y.subset <- Y[single.subs]
      if (anyNA(Y.subset)){
        stop("Estimate: Missing values in data.")
      }
    }
    if (!all(deterministic.list.newdata$is.deterministic |
             deterministic.g.list.newdata$is.deterministic)) {
      if (is.null(fit.and.predict) || multiple.Qstar || multiple.subs) {
        ## if(names(data)[cur.node] == "Y_3") {
        ## browser()
        ## print(names(data)[cur.node])
        ## }
        if(inputs$verbose)message("Regressing ",names(data)[cur.node]," on history with ",NROW(Y.subset)," observations")
        fit.and.predict <- FitAndPredict()
        m <- fit.and.predict$m
        predicted.values[, regime.index] <- fit.and.predict$predicted.values
      }
      else {
        predicted.values[, regime.index] <- PredictOnly(newdata)
      }
      ## print(calc.meanL)
      if (calc.meanL)
        prob.A.is.1.meanL[, regime.index, ] <- PredictProbAMeanL()
    }
    else {
      m <- "all rows are deterministic, no estimation took place"
    }
    predicted.values[deterministic.g.list.newdata$is.deterministic,
                     regime.index] <- deterministic.g.list.newdata$prob1
    if (calc.meanL)
      prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic,
                        regime.index, ] <- deterministic.g.list.newdata$prob1
    is.deterministic[, regime.index] <- deterministic.list.newdata$is.deterministic
    if (!called.from.estimate.g)
      deterministic.Q[deterministic.list.newdata$is.deterministic,
                      regime.index] <- deterministic.list.newdata$Q
    if (isTRUE(attr(SL.library, "return.fit", exact = TRUE))) {
      fit[[regime.index]] <- m
    }
    else {
      if (use.glm) {
        if (class(m)[1] %in% c("speedglm", "glm")) {
          fit[[regime.index]] <- summary(m)$coefficients
        }
        else {
          stopifnot(class(m)[1] %in% c("no.Y.variation",
                                       "character"))
          fit[[regime.index]] <- m
        }
      }
      else {
        if (inherits(m,"ltmle.glmnet")){
          fit[[regime.index]] <- m$selected_beta
        } else{
          capture.output(print.m <- print(m))
          fit[[regime.index]] <- print.m
        }
      }
    }
  }
  return(list(predicted.values = predicted.values, fit = fit,
              is.deterministic = is.deterministic, deterministic.Q = deterministic.Q,
              prob.A.is.1.meanL = prob.A.is.1.meanL))
}
EstimateG <- function (inputs)
{
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  nodes <- inputs$all.nodes
  g <- cum.g <- cum.g.unbounded <- prob.A.is.1 <- array(NaN,
                                                        dim = c(n, length(nodes$AC), num.regimes))
  if (inputs$variance.method == "ic") {
    cum.g.meanL <- cum.g.meanL.unbounded <- NULL
  }
  else {
    g.meanL <- cum.g.meanL <- cum.g.meanL.unbounded <- array(NaN,
                                                             dim = c(n, length(nodes$AC), num.regimes, length(nodes$LY) -
                                                                       1))
  }
  fit <- vector("list", length(nodes$AC))
  names(fit) <- names(inputs$data)[nodes$AC]
  if (inputs$variance.method != "ic" && anyNA(inputs$regimes)) {
    regimes.meanL <- inputs$regimes
    for (i in seq_along(nodes$A)) {
      for (regime.index in 1:num.regimes) {
        regimes.meanL[is.na(regimes.meanL[, i, regime.index]),
                      i, regime.index] <- Mode(inputs$regimes[, i,
                                                              regime.index], na.rm = TRUE)
      }
    }
  }
  else {
    regimes.meanL <- NULL
  }
  for (i in seq_along(nodes$AC)) {
    cur.node <- nodes$AC[i]
    uncensored <- IsUncensored(inputs$uncensored, nodes$C,
                               cur.node)
    deterministic.origdata <- IsDeterministic(inputs$data,
                                              cur.node, inputs$deterministic.Q.function, nodes,
                                              called.from.estimate.g = TRUE, inputs$survivalOutcome)$is.deterministic
    if (is.numeric(inputs$gform)) {
      if (!is.null(inputs$deterministic.g.function))
        stop("deterministic.g.function is not compatible with numeric gform")
      prob.A.is.1[, i, ] <- inputs$gform[, i, ]
      g.est <- list(is.deterministic = deterministic.origdata)
      fit[[i]] <- "no fit due to numeric gform"
    }
    else {
      form <- inputs$gform[i]
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data,
                                                        cur.node, inputs$deterministic.g.function, nodes,
                                                        using.newdata = F)
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic
      #
      # subsetting data at the current node
      #
      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs$intervention.match,
                                                nodes$A, cur.node = nodes$AC[i])
        subs <- uncensored & intervention.match & !deterministic.origdata &
          !deterministic.g.origdata
      }
      else {
        subs <- uncensored & !deterministic.origdata &
          !deterministic.g.origdata
      }
      if (inputs$verbose){ message("EstimateG: estimating formula ",form," ...")}
      g.est <- Estimate(inputs, form = form, Qstar.kplus1 = NULL,
                        subs = subs, family = quasibinomial(), type = "response",
                        nodes = nodes, called.from.estimate.g = TRUE,
                        calc.meanL = inputs$variance.method != "ic",
                        cur.node = cur.node, regimes.meanL = regimes.meanL,
                        regimes.with.positive.weight = 1:num.regimes)
      prob.A.is.1[, i, ] <- g.est$predicted.values
      fit[[i]] <- g.est$fit
    }
    if (cur.node %in% nodes$A) {
      cur.abar <- AsMatrix(inputs$regimes[, nodes$A ==
                                            cur.node, ])
      if (is.null(regimes.meanL)) {
        cur.abar.meanL <- cur.abar
      }
      else {
        cur.abar.meanL <- AsMatrix(regimes.meanL[, nodes$A ==
                                                   cur.node, ])
      }
    }
    else {
      cur.abar <- cur.abar.meanL <- matrix(1, nrow(inputs$data),
                                           num.regimes)
    }
    g[, i, ] <- CalcG(AsMatrix(prob.A.is.1[, i, ]), cur.abar,
                      g.est$is.deterministic)
    if (inputs$variance.method != "ic") {
      if (is.numeric(inputs$gform)) {
        if (anyNA(g[, i, ]))
          stop("Error - NA in numeric gform. There may not be NA values in gform (including after censoring if variance.method is 'tmle' or 'iptw'.")
        g.meanL[, i, , ] <- g[, i, ]
      }
      else {
        for (j in sseq(1, dim(g.meanL)[4])) {
          g.meanL[, i, , j] <- CalcG(AsMatrix(g.est$prob.A.is.1.meanL[,
                                                                      , j]), cur.abar.meanL, g.est$is.deterministic)
        }
      }
    }
    if (anyNA(g[uncensored, i, ]))
      stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")
  }
  for (regime.index in 1:num.regimes) {
    cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]),
                           inputs$gbounds)
    cum.g[, , regime.index] <- cum.g.list$bounded
    cum.g.unbounded[, , regime.index] <- cum.g.list$unbounded
    if (inputs$variance.method != "ic") {
      for (j in sseq(1, dim(g.meanL)[4])) {
        cum.g.list <- CalcCumG(AsMatrix(g.meanL[, , regime.index,
                                                j]), inputs$gbounds)
        cum.g.meanL[, , regime.index, j] <- cum.g.list$bounded
        cum.g.meanL.unbounded[, , regime.index, j] <- cum.g.list$unbounded
      }
    }
  }
  return(list(cum.g = cum.g, cum.g.unbounded = cum.g.unbounded,
              cum.g.meanL = cum.g.meanL, fit = ReorderFits(fit), prob.A.is.1 = prob.A.is.1,
              cum.g.meanL.unbounded = cum.g.meanL.unbounded))
}
EstimateTime <-
  function (inputs) 
  {
    sample.size <- 50
    if (nrow(inputs$data) < sample.size) {
      message(paste("Timing estimate unavailable when n <", 
                    sample.size))
      return(NULL)
    }
    sample.index <- sample(nrow(inputs$data), size = sample.size)
    small.inputs <- inputs
    small.inputs$data <- small.inputs$data[sample.index, ]
    small.inputs$data <- droplevels(small.inputs$data)
    small.inputs$regimes <- small.inputs$regimes[sample.index, 
                                                 , , drop = F]
    small.inputs$observation.weights <- small.inputs$observation.weights[sample.index]
    small.inputs$uncensored <- small.inputs$uncensored[sample.index, 
                                                       , drop = F]
    small.inputs$intervention.match <- small.inputs$intervention.match[sample.index, 
                                                                       , , drop = F]
    small.inputs$combined.summary.measures <- small.inputs$combined.summary.measures[sample.index, 
                                                                                     , , , drop = F]
    if (!is.null(small.inputs$id)) 
      small.inputs$id <- small.inputs$id[sample.index]
    if (is.numeric(inputs$gform)) 
      small.inputs$gform <- small.inputs$gform[sample.index, 
                                               , , drop = F]
    if (length(dim(inputs$msm.weights)) == 3) 
      small.inputs$msm.weights <- small.inputs$msm.weights[sample.index, 
                                                           , , drop = F]
    start.time <- Sys.time()
    try.result <- suppressWarnings(try(MainCalcs(small.inputs), 
                                       silent = TRUE))
    if (inherits(try.result, "try-error")) {
      message("Timing estimate unavailable")
    }
    else {
      elapsed.time <- Sys.time() - start.time
      est.time1 <- round(sqrt(as.double(elapsed.time, units = "mins") * 
                                nrow(inputs$data)/sample.size), digits = 0)
      est.time2 <- round(as.double(elapsed.time, units = "mins") * 
                           nrow(inputs$data)/sample.size, digits = 0)
      if (est.time2 == 0) {
        est.time.str <- "< 1 minute"
      }
      else if (est.time2 == 1) {
        est.time.str <- "1 minute"
      }
      else {
        est.time.str <- paste(est.time1, "to", est.time2, 
                              "minutes")
      }
      message("Estimate of time to completion: ", est.time.str)
    }
    return(NULL)
  }
EstimateVariance <-
  function (inputs, nodes, combined.summary.measures, regimes.with.positive.weight,
            uncensored, alive, Qstar, Qstar.kplus1, cur.node, msm.weights,
            LYnode.index, ACnode.index, cum.g, prob.A.is.1, cum.g.meanL,
            cum.g.unbounded, cum.g.meanL.unbounded, observation.weights,
            is.last.LYnode, intervention.match)
  {
    if (inputs$variance.method == "ic")
      return(NA)
    est.var.iptw <- inputs$variance.method == "iptw"
    TmleOfVariance <- function(Z, Z.meanL) {
      if (all(is.na(Z)))
        stop("all Z are NA in EstimateVariance")
      if (length(Z.meanL) == 0 || all(Z == 0 | is.na(Z))) {
        Qstar <- Scale(Z, 0, 1)
        return(list(EZd1 = mean(Z, na.rm = T), Qstar = Qstar))
      }
      if (est.var.iptw) {
        index <- uncensored & intervention.match[, d1]
        g <- cum.g[index, ACnode.index, d1]
        Y <- Scale(Z, 0, 1)[index]
        iptw.estimate <- sum(Y/g)/sum(1/g)
        return(list(EZd1 = iptw.estimate * diff(range(Z,
                                                      na.rm = T)) + min(Z, na.rm = T), Qstar = rep(NA,
                                                                                                   length(Z))))
      }
      sparsity.data <- inputs$data[, 1:cur.node]
      sparsity.data[, cur.node] <- Scale(Z, 0, 1)
      temp.nodes <- lapply(nodes, function(x) x[x <= cur.node])
      if (cur.node %in% temp.nodes$L) {
        temp.nodes$L <- setdiff(temp.nodes$L, cur.node)
        temp.nodes$Y <- c(temp.nodes$Y, cur.node)
      }
      stratify <- FALSE
      Qform <- paste(GetDefaultForm(sparsity.data[, 1:cur.node],
                                    nodes = temp.nodes, is.Qform = TRUE, stratify = stratify,
                                    survivalOutcome = FALSE, showMessage = FALSE), paste0("+ sparityAdj_Z.meanL_",
                                                                                          1:length(temp.nodes$LY)))
      Qform[length(Qform)] <- "IDENTITY"
      Z.meanL <- apply(AsMatrix(Z.meanL), 2, LogitScale)
      sparsity.data <- cbind(Z.meanL, sparsity.data)
      names(sparsity.data)[sseq(1, ncol(Z.meanL))] <- paste0("sparityAdj_Z.meanL_",
                                                             sseq(1, ncol(Z.meanL)))
      temp.nodes <- lapply(temp.nodes, function(x) x + ncol(Z.meanL))
      names(Qform) <- names(sparsity.data)[temp.nodes$LY]
      attr(sparsity.data, "called.from.estimate.variance") <- TRUE
      var.tmle <- Ltmle(sparsity.data, Anodes = temp.nodes$A,
                        Cnodes = temp.nodes$C, Lnodes = temp.nodes$L, Ynodes = temp.nodes$Y,
                        survivalOutcome = FALSE, Qform = Qform, gform = drop3(prob.A.is.1[,
                                                                                          1:ACnode.index, d1, drop = FALSE]), abar = GetABar(inputs$regimes,
                                                                                                                                             d1, temp.nodes$A), gbounds = inputs$gbounds,
                        stratify = stratify, estimate.time = FALSE, deterministic.Q.function = det.q.function,
                        variance.method = "ic", observation.weights = observation.weights)
      EZd1 <- var.tmle$estimates["tmle"] * diff(range(Z, na.rm = T)) +
        min(Z, na.rm = T)
      return(list(EZd1 = EZd1, Qstar = var.tmle$Qstar))
    }
    EqualRegimesIndex <- function(dd1, dd2) {
      if (!any(nodes$A <= cur.node))
        return(rep(TRUE, n))
      eq <- rowAlls(AsMatrix(inputs$regimes[, which(nodes$A <=
                                                      cur.node), dd1]) == AsMatrix(inputs$regimes[, which(nodes$A <=
                                                                                                            cur.node), dd2]))
      eq[is.na(eq)] <- FALSE
      return(eq)
    }
    IsStaticTreatment <- function() {
      for (dd1 in regimes.with.positive.weight) {
        for (dd2 in regimes.with.positive.weight[regimes.with.positive.weight >
                                                 dd1]) {
          if (any(EqualRegimesIndex(dd1, dd2)))
            return(FALSE)
        }
      }
      return(TRUE)
    }
    num.regimes <- dim(inputs$regimes)[3]
    num.betas <- ncol(combined.summary.measures)
    n <- nrow(inputs$data)
    if (inputs$survivalOutcome) {
      det.q.function <- function(data, current.node, nodes,
                                 called.from.estimate.g) {
        if (!any(nodes$Y < current.node))
          return(NULL)
        prev.Y <- data[, nodes$Y[nodes$Y < current.node],
                       drop = F]
        prev.Y[is.na(prev.Y)] <- 0
        is.deterministic <- rowAnys(prev.Y == 1)
        Q.value <- data[is.deterministic, max(nodes$Y)]
        return(list(is.deterministic = is.deterministic,
                    Q.value = Q.value))
      }
    }
    else {
      det.q.function <- NULL
    }
    static.treatment <- IsStaticTreatment()
    variance.estimate <- matrix(0, num.betas, num.betas)
    Sigma <- array(dim = c(n, num.regimes, num.regimes))
    if (!is.last.LYnode)
      Q.data <- inputs$data[alive, 1:cur.node, drop = F]
    for (d1 in regimes.with.positive.weight) {
      if (static.treatment) {
        d2.regimes <- d1
      }
      else {
        d2.regimes <- regimes.with.positive.weight
      }
      for (d2 in d2.regimes) {
        if (is.last.LYnode) {
          Sigma[, d1, d2] <- Qstar[, d1] * (1 - Qstar[,
                                                      d1])
        }
        else {
          if (any(alive)) {
            resid.sq <- (Qstar.kplus1[alive, d1] - Qstar[alive,
                                                         d1]) * (Qstar.kplus1[alive, d2] - Qstar[alive,
                                                                                                 d2])
            resid.sq.range <- range(resid.sq, na.rm = T)
            if (diff(resid.sq.range) > 0.0001) {
              Q.data[, cur.node] <- (resid.sq - resid.sq.range[1])/diff(resid.sq.range)
              names(Q.data)[cur.node] <- "Q.kplus1"
              m <- ltmle.glm(formula = formula(inputs$Qform[LYnode.index]),
                             family = quasibinomial(), data = Q.data,
                             weights = NULL)
              Q.newdata <- SetA(data = Q.data, regimes = inputs$regimes[alive,
                                                                        , d1, drop = F], Anodes = nodes$A, cur.node = cur.node)
              SuppressGivenWarnings(Q.resid.sq.pred <- predict(m,
                                                               newdata = Q.newdata, type = "response"),
                                    "prediction from a rank-deficient fit may be misleading")
              Sigma[alive, d1, d2] <- Q.resid.sq.pred *
                diff(resid.sq.range) + resid.sq.range[1]
            }
            else {
              resid.sq.value <- min(resid.sq, na.rm = T)
              Sigma[alive, d1, d2] <- resid.sq.value
            }
          }
          Sigma[!alive, d1, d2] <- 0
        }
      }
    }
    if (est.var.iptw)
      Z.without.sum.meas.meanL <- Z.meanL <- NA
    no.V <- length(inputs$baseline.column.names) == 0
    if ((!est.var.iptw && static.treatment) || (est.var.iptw &&
                                                static.treatment && no.V)) {
      for (d1 in regimes.with.positive.weight) {
        Z.without.sum.meas <- Sigma[, d1, d1]/cum.g[, ACnode.index,
                                                    d1] * cum.g.unbounded[, ACnode.index, d1]/cum.g[,
                                                                                                    ACnode.index, d1] * msm.weights[, d1]^2 * observation.weights^2
        if (!est.var.iptw)
          Z.without.sum.meas.meanL <- 1/cum.g.meanL[, ACnode.index,
                                                    d1, ] * cum.g.meanL.unbounded[, ACnode.index,
                                                                                  d1, ]/cum.g.meanL[, ACnode.index, d1, ] * msm.weights[,
                                                                                                                                        d1]^2 * observation.weights^2
        var.tmle <- TmleOfVariance(Z.without.sum.meas, Z.without.sum.meas.meanL)
        if (no.V) {
          variance.estimate <- variance.estimate + (combined.summary.measures[1,
                                                                              , d1] %*% t(combined.summary.measures[1, ,
                                                                                                                    d1])) * var.tmle$EZd1
        }
        else {
          baseline.msm <- paste("Qstar ~", paste(inputs$baseline.column.names,
                                                 collapse = " + "), "+", paste0("I(", inputs$baseline.column.names,
                                                                                "^2)", collapse = " + "))
          V.data <- data.frame(Qstar = var.tmle$Qstar,
                               inputs$data[, inputs$baseline.column.names,
                                           drop = FALSE])
          m <- ltmle.glm(formula(baseline.msm), family = quasibinomial(),
                         data = V.data, weights = NULL)
          SuppressGivenWarnings(pred.Qstar <- predict(m,
                                                      type = "response", newdata = V.data) * diff(range(Z.without.sum.meas,
                                                                                                        na.rm = T)) + min(Z.without.sum.meas, na.rm = T),
                                "prediction from a rank-deficient fit may be misleading")
          variance.estimate.sum <- crossprod(combined.summary.measures[,
                                                                       , d1], combined.summary.measures[, , d1] *
                                               pred.Qstar)
          variance.estimate <- variance.estimate + variance.estimate.sum/n
        }
      }
    }
    else {
      for (beta.index2 in 1:num.betas) {
        for (d1 in regimes.with.positive.weight) {
          Z.base <- rep(0, n)
          if (!est.var.iptw)
            Z.base.meanL <- matrix(0, n, dim(cum.g.meanL)[4])
          for (d2 in regimes.with.positive.weight) {
            equal.regimes.index <- EqualRegimesIndex(d1,
                                                     d2)
            h1 <- combined.summary.measures[, beta.index2,
                                            d2] * msm.weights[, d2]
            Z.base[equal.regimes.index] <- Z.base[equal.regimes.index] +
              h1[equal.regimes.index] * Sigma[equal.regimes.index,
                                              d1, d2]/cum.g[equal.regimes.index, ACnode.index,
                                                            d1] * observation.weights[equal.regimes.index]
            if (!est.var.iptw)
              Z.base.meanL[equal.regimes.index, ] <- Z.base.meanL[equal.regimes.index,
              ] + h1[equal.regimes.index] * 1/cum.g.meanL[equal.regimes.index,
                                                          ACnode.index, d1, ] * observation.weights[equal.regimes.index]
          }
          for (beta.index1 in 1:num.betas) {
            if (beta.index1 >= beta.index2) {
              Z <- combined.summary.measures[, beta.index1,
                                             d1] * msm.weights[, d1] * cum.g.unbounded[,
                                                                                       ACnode.index, d1]/cum.g[, ACnode.index,
                                                                                                               d1] * observation.weights * Z.base
              if (!est.var.iptw)
                Z.meanL <- combined.summary.measures[,
                                                     beta.index1, d1] * msm.weights[, d1] *
                  cum.g.meanL.unbounded[, ACnode.index,
                                        d1, ]/cum.g.meanL[, ACnode.index, d1,
                                        ] * observation.weights * Z.base.meanL
              var.tmle <- TmleOfVariance(Z, Z.meanL)
              variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index1,
                                                                               beta.index2] + var.tmle$EZd1
            }
            else {
              variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index2,
                                                                               beta.index1]
            }
          }
        }
      }
    }
    if (max(abs(variance.estimate - t(variance.estimate))) >
        1.0000000000000001e-05)
      stop("not symmetric")
    if (any(eigen(variance.estimate, only.values = TRUE)$values <
            -1e-08)) {
      variance.estimate <- MakePosDef(variance.estimate)
    }
    return(variance.estimate)
  }
FinalizeIC <-
  function (IC, combined.summary.measures, Qstar, m.beta, msm.weights, 
            observation.weights, id) 
  {
    num.betas <- ncol(IC)
    n <- nrow(Qstar)
    num.regimes <- ncol(Qstar)
    num.final.Ynodes <- dim(Qstar)[3]
    stopifnot(num.betas == ncol(combined.summary.measures))
    finalIC <- matrix(0, nrow = n, ncol = num.betas)
    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimes) {
        if (any(msm.weights[, i, j] > 0)) {
          m1 <- matrix(Qstar[, i, j] - m.beta[, i, j], 
                       ncol = 1)
          for (k in 1:num.betas) {
            m2 <- combined.summary.measures[, k, i, j]
            finalIC[, k] <- finalIC[, k] + msm.weights[, 
                                                       i, j] * observation.weights * (m1 * m2)
          }
        }
      }
    }
    IC <- IC + finalIC
    return(HouseholdIC(IC, id))
  }
FitPooledMSM <-
  function (working.msm, Qstar, combined.summary.measures, msm.weights) 
  {
    n <- dim(Qstar)[1]
    num.regimes <- dim(Qstar)[2]
    num.final.Ynodes <- dim(Qstar)[3]
    num.summary.measures <- dim(combined.summary.measures)[2]
    X <- apply(combined.summary.measures, 2, rbind)
    Y <- as.vector(Qstar)
    weight.vec <- as.vector(msm.weights)
    data.pooled <- data.frame(Y, X)
    positive.weight <- weight.vec > 0
    m <- ltmle.glm(formula(working.msm), data = data.pooled[positive.weight, 
    ], family = quasibinomial(), weights = weight.vec[positive.weight])
    SuppressGivenWarnings(m.beta <- predict(m, newdata = data.pooled, 
                                            type = "response"), "prediction from a rank-deficient fit may be misleading")
    dim(m.beta) <- dim(Qstar)
    return(list(m = m, m.beta = m.beta))
  }
FixedTimeTMLE <-
  function (inputs, nodes, msm.weights, combined.summary.measures,
            g.list)
  {
    data <- inputs$data
    num.regimes <- dim(inputs$regimes)[3]
    n <- nrow(data)
    num.betas <- ncol(combined.summary.measures)
    tmle <- rep(NA, num.regimes)
    IC <- matrix(0, nrow = n, ncol = num.betas)
    cum.g.used <- array(FALSE, dim = dim(g.list$cum.g))
    est.var <- matrix(0, num.betas, num.betas)
    regimes.with.positive.weight <- which(apply(msm.weights >
                                                  0, 2, any))
    if (length(regimes.with.positive.weight) == 0)
      stop("All regimes have weight 0 (one possible reason is that msm.weights='emipirical' and no data rows match any of the regimes and are uncensored)")
    fit.Qstar <- fit.Q <- vector("list", length(nodes$LY))
    names(fit.Qstar) <- names(fit.Q) <- names(data)[nodes$LY]
    Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow = n, ncol = num.regimes)
    mean.summary.measures <- apply(abs(combined.summary.measures),
                                   2, mean)
    if (length(nodes$LY) > 0) {
      for (LYnode.index in length(nodes$LY):1) {
        cur.node <- nodes$LY[LYnode.index]
        deterministic.list.origdata <- IsDeterministic(data,
                                                       cur.node, inputs$deterministic.Q.function, nodes,
                                                       called.from.estimate.g = FALSE, inputs$survivalOutcome)
        uncensored <- IsUncensored(inputs$uncensored, nodes$C,
                                   cur.node)
        intervention.match <- InterventionMatch(inputs$intervention.match,
                                                nodes$A, cur.node)
        if (inputs$stratify) {
          subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic
        }
        else {
          subs <- uncensored & !deterministic.list.origdata$is.deterministic
        }
        ## print(inputs$Qform[LYnode.index])
        if (inputs$verbose){
          message("FixedTimeTMLE: estimating Q for ",
                  names(data)[cur.node],
                  " using ",
                  sum(subs),
                  " observations.")
          ## output = c(output,list(data.table(Outcome = names(data)[cur.node],Atrisk = sum(subs))))
        }
        Q.est <- Estimate(inputs, form = inputs$Qform[LYnode.index],
                          Qstar.kplus1 = if (LYnode.index == length(nodes$LY))
                            Qstar.kplus1[, 1]
                          else Qstar.kplus1, family = quasibinomial(),
                          subs = subs, type = "link", nodes = nodes, called.from.estimate.g = FALSE,
                          calc.meanL = FALSE, cur.node = cur.node, regimes.meanL = NULL,
                          regimes.with.positive.weight = regimes.with.positive.weight)
        logitQ <- Q.est$predicted.values
        fit.Q[[LYnode.index]] <- Q.est$fit
        ACnode.index <- which.max(nodes$AC[nodes$AC < cur.node])
        SuppressGivenWarnings(update.list <- UpdateQ(Qstar.kplus1,
                                                     logitQ, combined.summary.measures, g.list$cum.g[,
                                                                                                     ACnode.index, ], inputs$working.msm, uncensored,
                                                     intervention.match, deterministic.list.origdata$is.deterministic,
                                                     msm.weights, inputs$gcomp, inputs$observation.weights),
                              GetWarningsToSuppress(update.step = TRUE))
        if (length(ACnode.index) > 0)
          cum.g.used[, ACnode.index, ] <- cum.g.used[,
                                                     ACnode.index, ] | update.list$cum.g.used
        Qstar <- update.list$Qstar
        Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic]
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
                        uncensored, intervention.match, regimes.with.positive.weight)
        curIC.relative.error <- abs(colSums(curIC))
        curIC.relative.error[mean.summary.measures > 0] <- curIC.relative.error[mean.summary.measures >
                                                                                  0]/mean.summary.measures[mean.summary.measures >
                                                                                                             0]
        if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {
          SetSeedIfRegressionTesting()
          fix.score.list <- FixScoreEquation(Qstar.kplus1,
                                             update.list$h.g.ratio, uncensored, intervention.match,
                                             Q.est$is.deterministic, Q.est$deterministic.Q,
                                             update.list$off, update.list$X, regimes.with.positive.weight)
          Qstar <- fix.score.list$Qstar
          curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
                          uncensored, intervention.match, regimes.with.positive.weight)
          update.list$fit <- fix.score.list$fit
        }
        est.var <- est.var + EstimateVariance(inputs, nodes,
                                              combined.summary.measures, regimes.with.positive.weight,
                                              uncensored, alive = !deterministic.list.origdata$is.deterministic,
                                              Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index,
                                              ACnode.index, g.list$cum.g, g.list$prob.A.is.1,
                                              g.list$cum.g.meanL, g.list$cum.g.unbounded, g.list$cum.g.meanL.unbounded,
                                              inputs$observation.weights, is.last.LYnode = (LYnode.index ==
                                                                                              length(nodes$LY)), intervention.match)
        IC <- IC + curIC
        Qstar.kplus1 <- Qstar
        fit.Qstar[[LYnode.index]] <- update.list$fit
      }
    }
    else {
      Qstar <- Qstar.kplus1
    }
    return(list(IC = IC, Qstar = Qstar, est.var = est.var, fit = list(Q = ReorderFits(fit.Q),
                                                                      Qstar = fit.Qstar), cum.g.used = cum.g.used))
  }
FixScoreEquation <-
  function (Qstar.kplus1, h.g.ratio, uncensored, intervention.match, 
            is.deterministic, deterministic.Q, off, X, regimes.with.positive.weight) 
  {
    CalcScore <- function(e) {
      Qstar <- QstarFromE(e)
      ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, 
                       intervention.match, regimes.with.positive.weight)
      return(sum(colSums(ICtemp)^2))
    }
    QstarFromE <- function(e) {
      Qstar <- plogis(off + X %*% e)
      dim(Qstar) <- dim(Qstar.kplus1)
      Qstar[is.deterministic] <- deterministic.Q[is.deterministic]
      return(Qstar)
    }
    FindMin <- function(minimizer) {
      num.tries <- 20
      init.e <- numeric(num.betas)
      for (i in 1:num.tries) {
        m <- nlminb(start = init.e, objective = CalcScore, 
                    control = list(abs.tol = max.objective, eval.max = 500, 
                                   iter.max = 500, x.tol = 1e-14, rel.tol = 1e-14))
        e <- m$par
        obj.val <- m$objective
        if (obj.val < max.objective) {
          m$ltmle.msg <- "updating step using glm failed to solve score equation; solved using nlminb"
          return(list(e = e, solved = TRUE, m = m))
        }
        init.e <- rnorm(num.betas)
      }
      return(list(e = numeric(num.betas), solved = FALSE, m = "score equation not solved!"))
    }
    max.objective <- 0.0001^2
    num.betas <- ncol(X)
    for (offset.lbound in c(1e-08, 0.0001, 0.001, 0.01)) {
      off <- Bound(off, qlogis(c(offset.lbound, 1 - offset.lbound)))
      l <- FindMin("nlminb")
      if (l$solved) 
        break
    }
    if (!l$solved) 
      stop("minimizer failed")
    Qstar <- QstarFromE(l$e)
    return(list(Qstar = Qstar, fit = l$m))
  }
GetABar <-
  function (regimes, regime.index, Anodes) 
  {
    abar <- AsMatrix(regimes[, seq_along(Anodes), regime.index])
    return(abar)
  }
GetCI <-
  function (estimate, std.dev, n) 
  {
    if (n < 100) {
      x <- qt(0.97499999999999998, df = n - 1) * std.dev
    }
    else {
      x <- qnorm(0.97499999999999998) * std.dev
    }
    CI <- cbind(`2.5%` = estimate - x, `97.5%` = estimate + x)
    return(CI)
  }
GetDefaultForm <-
  function (data, nodes, is.Qform, stratify, survivalOutcome, showMessage) 
  {
    if (is.Qform) {
      lhs <- rep("Q.kplus1", length(nodes$LY))
      node.set <- nodes$LY
    }
    else {
      lhs <- names(data)[nodes$AC]
      node.set <- nodes$AC
    }
    if (stratify) {
      stratify.nodes <- c(nodes$C, nodes$A)
    }
    else {
      stratify.nodes <- c(nodes$C)
    }
    if (survivalOutcome) {
      stratify.nodes <- c(stratify.nodes, nodes$Y)
    }
    form <- NULL
    for (i in seq_along(node.set)) {
      cur.node <- node.set[i]
      if (cur.node == 1) {
        form[i] <- paste(lhs[i], "~ 1")
      }
      else {
        parent.node.names <- names(data)[setdiff(1:(cur.node - 
                                                      1), stratify.nodes)]
        if (length(parent.node.names) == 0) {
          form[i] <- paste(lhs[i], "~ 1")
        }
        else {
          form[i] <- paste(lhs[i], "~", paste(parent.node.names, 
                                              collapse = " + "))
        }
      }
      names(form)[i] <- names(data)[cur.node]
    }
    if (showMessage) {
      message(ifelse(is.Qform, "Qform", "gform"), " not specified, using defaults:")
      lapply(seq_along(form), function(i, names) {
        message("formula for ", names[i], ":")
        message(capture.output(print(as.formula(form[i]), 
                                     showEnv = FALSE)))
      }, names = names(form))
      message("")
    }
    return(form)
  }
GetLibrary <-
  function (SL.library, estimate.type) 
  {
    NullToGlm <- function(libr) if (is.null(libr)) 
      "glm"
    else libr
    if (is.null(names(SL.library))) 
      return(NullToGlm(SL.library))
    if (!identical(sort(names(SL.library)), sort(c("Q", "g")))) 
      stop("If SL.library has names, it must have two names: Q and g")
    if (!estimate.type %in% c("Q", "g")) 
      stop("bad estimate.type")
    if (length(setdiff(names(attributes(SL.library)), c("names", 
                                                        "return.fit"))) > 0) 
      stop("If SL.library has attributes, the only valid attributes are name and return.fit")
    lib <- SL.library[[estimate.type]]
    attr(lib, "return.fit") <- attr(SL.library, "return.fit", 
                                    exact = TRUE)
    return(NullToGlm(SL.library[[estimate.type]]))
  }
GetMSMInputsForLtmle <-
  function (data, abar, rule, gform) 
  {
    if ((!missing(abar) && is.list(abar)) || is.list(rule)) {
      if (is.list(rule)) {
        if (length(rule) != 2) 
          stop("If rule is a list, it must be of length 2")
        regimes1 <- RegimesFromAbar(data, abar, rule[[1]])
        regimes0 <- RegimesFromAbar(data, abar, rule[[2]])
      }
      else {
        if (length(abar) != 2) 
          stop("If abar is a list, it must be of length 2")
        regimes1 <- RegimesFromAbar(data, abar[[1]], rule)
        regimes0 <- RegimesFromAbar(data, abar[[2]], rule)
      }
      if (ncol(regimes1) != ncol(regimes0)) 
        stop("If abar or rule is a list, both elements must give a matrix with the same number of columns")
      if (nrow(regimes1) != nrow(regimes0)) 
        stop("If abar or rule is a list, both elements must give a matrix with the same number of rows")
      regimes <- c(regimes1, regimes0)
      dim(regimes) <- c(nrow(regimes1), ncol(regimes1), 2)
      summary.measures <- array(1:0, dim = c(2, 1, 1))
      colnames(summary.measures) <- "A"
      working.msm <- "Y ~ A"
      msm.weights <- matrix(1, nrow = 2, ncol = 1)
    }
    else {
      regimes <- RegimesFromAbar(data, abar, rule)
      working.msm <- "Y ~ 1"
      msm.weights <- matrix(1, nrow = 1, ncol = 1)
      summary.measures <- array(dim = c(1, 0, 1))
    }
    msm.inputs <- list(regimes = regimes, working.msm = working.msm, 
                       summary.measures = summary.measures, gform = gform, final.Ynodes = NULL, 
                       msm.weights = msm.weights)
    return(msm.inputs)
  }
GetMsmWeights <-
  function (inputs) 
  {
    n <- nrow(inputs$data)
    num.regimes <- dim(inputs$regimes)[3]
    stopifnot(num.regimes >= 1)
    num.final.Ynodes <- length(inputs$final.Ynodes)
    if (is.equal(inputs$msm.weights, "empirical")) {
      msm.weights <- matrix(nrow = num.regimes, ncol = num.final.Ynodes)
      for (j in 1:num.final.Ynodes) {
        final.Ynode <- inputs$final.Ynodes[j]
        regimes.subset <- inputs$regimes[, inputs$all.nodes$A < 
                                           final.Ynode, , drop = FALSE]
        if (ncol(regimes.subset) > 0) {
          is.duplicate <- duplicated(regimes.subset, MARGIN = 3)
        }
        else {
          is.duplicate <- c(FALSE, rep(TRUE, num.regimes - 
                                         1))
        }
        uncensored <- IsUncensored(inputs$uncensored, inputs$all.nodes$C, 
                                   cur.node = final.Ynode)
        intervention.match <- InterventionMatch(inputs$intervention.match, 
                                                inputs$all.nodes$A, cur.node = final.Ynode)
        for (i in 1:num.regimes) {
          if (is.duplicate[i]) {
            msm.weights[i, j] <- 0
          }
          else {
            msm.weights[i, j] <- sum(uncensored & intervention.match[, 
                                                                     i])/nrow(inputs$data)
          }
        }
      }
    }
    else if (is.null(inputs$msm.weights)) {
      msm.weights <- array(1, dim = c(n, num.regimes, num.final.Ynodes))
    }
    else {
      msm.weights <- inputs$msm.weights
    }
    if (is.equal(dim(msm.weights), c(num.regimes, num.final.Ynodes))) {
      msm.weights <- array(rep(msm.weights, each = n), dim = c(n, 
                                                               num.regimes, num.final.Ynodes))
    }
    if (anyNA(msm.weights) || any(msm.weights < 0)) 
      stop("all msm.weights must be >= 0 and not NA")
    return(msm.weights)
  }
GetPValue <-
  function (q, n) 
  {
    if (n < 100) {
      return(pt(q, n - 1))
    }
    else {
      return(pnorm(q))
    }
  }
GetSummary <-
  function (eff.list, cov.mat, n) 
  {
    estimate <- eff.list$est
    v <- t(eff.list$gradient) %*% cov.mat %*% eff.list$gradient
    stopifnot(length(v) == 1)
    std.dev <- sqrt(v[1, 1]/n)
    if (eff.list$log.std.err) {
      pvalue <- 2 * GetPValue(-abs(log(estimate)/std.dev), 
                              n)
      CI <- exp(GetCI(log(estimate), std.dev, n))
    }
    else {
      pvalue <- 2 * GetPValue(-abs(estimate/std.dev), n)
      CI <- GetCI(estimate, std.dev, n)
    }
    CI <- Bound(CI, eff.list$CIBounds)
    return(list(long.name = eff.list$long.name, estimate = estimate, 
                std.dev = std.dev, pvalue = pvalue, CI = CI, log.std.err = eff.list$log.std.err))
  }
GetSummaryLtmleMSMInfo <-
  function (object, estimator) 
  {
    if (!estimator %in% c("tmle", "iptw", "gcomp")) 
      stop("estimator should be one of: tmle, iptw, gcomp")
    if (estimator == "tmle") {
      if (object$gcomp) 
        stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
      estimate <- object$beta
      IC <- object$IC
    }
    else if (estimator == "iptw") {
      estimate <- object$beta.iptw
      IC <- object$IC.iptw
    }
    else if (estimator == "gcomp") {
      if (!object$gcomp) 
        stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
      estimate <- object$beta
      IC <- object$IC
    }
    IC.variance <- apply(IC, 2, var)
    if (is.null(object$variance.estimate)) {
      v <- IC.variance
    }
    else {
      v <- pmax(diag(object$variance.estimate), IC.variance)
    }
    variance.estimate.ratio <- v/IC.variance
    return(list(estimate = estimate, IC = IC, variance.estimate.ratio = variance.estimate.ratio, 
                v = v))
  }
GetWarningsToSuppress <-
  function (update.step = FALSE) 
  {
    warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                              "prediction from a rank-deficient fit may be misleading")
    if (update.step) {
      warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
    }
    return(warnings.to.suppress)
  }
HouseholdIC <-
  function (recordIC, id) 
  {
    if (is.null(id)) 
      return(recordIC)
    householdIC <- as.matrix(aggregate(recordIC, list(id = id), 
                                       sum)[, -1, drop = FALSE])
    num.records <- nrow(recordIC)
    num.households <- nrow(householdIC)
    householdIC <- householdIC * num.households/num.records
    return(householdIC)
  }
InterventionMatch <-
  function (intervention.match.array, Anodes, cur.node) 
  {
    index <- which.max(Anodes[Anodes < cur.node])
    if (length(index) == 0) 
      return(matrix(TRUE, nrow(intervention.match.array), dim(intervention.match.array)[3]))
    return(AsMatrix(intervention.match.array[, index, ]))
  }
is.equal <-
  function (...) 
  {
    isTRUE(all.equal(...))
  }
is.glm <-
  function (SL.library) 
  {
    is.equal(SL.library, "glm", check.attributes = FALSE)
  }
IsBinary <-
  function (mat) 
  {
    is.equal(mat, as.numeric(as.logical(mat)))
  }
IsDeterministic <-
  function (data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, 
            survivalOutcome) 
  {
    if (survivalOutcome && any(nodes$Y < cur.node)) {
      last.Ynode <- max(nodes$Y[nodes$Y < cur.node])
      is.deterministic <- data[, last.Ynode] %in% TRUE
    }
    else {
      is.deterministic <- rep(FALSE, nrow(data))
    }
    default <- list(is.deterministic = is.deterministic, Q.value = 1)
    if (is.null(deterministic.Q.function)) 
      return(default)
    det.list <- deterministic.Q.function(data = data, current.node = cur.node, 
                                         nodes = nodes, called.from.estimate.g = called.from.estimate.g)
    if (is.null(det.list)) 
      return(default)
    if (called.from.estimate.g) {
      if (!is.list(det.list) || !("is.deterministic" %in% names(det.list)) || 
          !(length(det.list) %in% 1:2)) 
        stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
    }
    else {
      if (!is.list(det.list) || !setequal(names(det.list), 
                                          c("is.deterministic", "Q.value")) || length(det.list) != 
          2) 
        stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
    }
    if (!length(det.list$Q.value) %in% c(1, length(which(det.list$is.deterministic)))) 
      stop("the length of the 'Q.value' element of deterministic.Q.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
    Q.value.from.function <- rep(NA, nrow(data))
    Q.value.from.function[det.list$is.deterministic] <- det.list$Q.value
    set.by.function.and.death <- is.deterministic & det.list$is.deterministic
    if (any(Q.value.from.function[set.by.function.and.death] != 
            1)) {
      stop(paste("inconsistent deterministic Q at node:", names(data)[cur.node]))
    }
    finalY <- data[, max(nodes$Y)]
    # FIXME: the original code compared a singleton against a vector
    #browser()
    inconsistent.rows <- (det.list$Q.value %in% c(0, 1)) &
      (det.list$Q.value != finalY[det.list$is.deterministic]) &
      !is.na(finalY[det.list$is.deterministic])
    if (any(inconsistent.rows)){
      stop(paste("At node:", names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", 
                 paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), 
                       collapse = " ")))
    }
    Q.value <- rep(NA, nrow(data))
    Q.value[is.deterministic] <- 1
    Q.value[det.list$is.deterministic] <- det.list$Q.value
    is.deterministic <- is.deterministic | det.list$is.deterministic
    Q.value <- Q.value[is.deterministic]
    if (anyNA(is.deterministic) || anyNA(Q.value)) 
      stop("NA in is.deterministic or Q.value")
    return(list(is.deterministic = is.deterministic, Q.value = Q.value))
  }
IsDeterministicG <-
  function (data, cur.node, deterministic.g.function, nodes, using.newdata) 
  {
    stopifnot(cur.node %in% nodes$AC)
    default <- list(is.deterministic = rep(FALSE, nrow(data)), 
                    prob1 = NULL)
    if (is.null(deterministic.g.function)) 
      return(default)
    det.list <- deterministic.g.function(data = data, current.node = cur.node, 
                                         nodes = nodes)
    if (is.null(det.list)) 
      return(default)
    if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", 
                                                           "prob1")) || length(det.list) != 2) 
      stop("deterministic.g.function should return a list with names: is.deterministic, prob1")
    if (!length(det.list$prob1) %in% c(1, length(which(det.list$is.deterministic)))) 
      stop("the length of the 'prob1' element of deterministic.g.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
    inconsistent.rows <- (det.list$prob1 %in% c(0, 1)) & (det.list$prob1 != 
                                                            data[det.list$is.deterministic, cur.node]) & !is.na(data[det.list$is.deterministic, 
                                                                                                                     cur.node])
    if (any(inconsistent.rows)) {
      err.msg <- paste("At node:", names(data)[cur.node], "deterministic.g.function is inconsistent with data - prob1 is either 0 or 1 but this does not match the node value.\nCheck data rows:", 
                       paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), 
                             collapse = " "))
      if (using.newdata) {
        err.msg <- paste(err.msg, "\n This error occured while calling deterministic.g.function on data where Anodes are set to abar.")
        cat("deterministic.g.function is inconsistent with data.\nAfter setting Anodes to abar, the data looks like this:\n")
        print(head(data[det.list$is.deterministic[inconsistent.rows], 
        ]))
      }
      stop(err.msg)
    }
    return(det.list)
  }
IsUncensored <-
  function (uncensored.matrix, Cnodes, cur.node) 
  {
    index <- which.max(Cnodes[Cnodes < cur.node])
    if (length(index) == 0) 
      return(rep(TRUE, nrow(uncensored.matrix)))
    return(uncensored.matrix[, index])
  }
LhsVars <-
  function (f) 
  {
    f <- as.formula(f)
    return(as.character(f[[2]]))
  }
LogitScale <-
  function (x) 
  {
    qlogis(Scale(x, 0.01, 0.98999999999999999))
  }
#
# Function called from Estimate (via EstimateG and FixedTimeTMLE)
#
ltmle.glm.fit <- function (x, y, weights, family, offset, intercept) {
  ## print(table(y))
  ## print(colnames(y))
  if (all(weights==1)||is.null(weights)){
    try.speedglm <- try({
      m <- speedglm::speedglm.wfit(y = y,
                                   X = x,
                                   family = family,
                                   offset = offset,
                                   intercept = intercept,
                                   maxit = 100)
      class(m) <- c("speedglm", "speedlm")
    }, silent = TRUE)
    if (inherits(try.speedglm, "try-error")) {
      ShowGlmMessage()
      try.glm <- try({m <- glm.fit(x = x,
                                   y = y,
                                   family = family,
                                   offset = offset,
                                   intercept = intercept,
                                   control = glm.control(maxit = 100))
      })
      if (inherits(try.glm,"try-error")){
        stop("Could not fit glm")
      }
      class(m) <- c("glm", "lm")
    }
  }else{
    try.glm <- try({
      m <- glm.fit(x = x,
                   y = y,
                   family = family,
                   weights = weights,
                   offset = offset,
                   intercept = intercept,
                   control = glm.control(maxit = 100))})
    if (inherits(try.glm,"try-error")) stop("Could not fit glm")
    class(m) <- c("glm", "lm")
  }
  return(m)
}

ltmle.sg <-
  function (d, Inodes, Lnodes, Ynodes, Qform, gform, gbd = 0, family = "quasibinomial", 
            move.to.weight) 
  {
    estQ <- function(Q.kplus1, d, Qform, uncensored, deterministic, 
                     h = NULL, family) {
      Qform <- update.formula(Qform, Q.kplus1 ~ .)
      m <- glm(as.formula(Qform), data = data.frame(Q.kplus1, 
                                                    d)[uncensored & !deterministic, ], family = family)
      Q1W <- predict(m, newdata = d, type = "response")
      if (!is.null(h)) {
        off <- qlogis(Bound(Q1W, c(0.0001, 0.99990000000000001)))
        if (move.to.weight) {
          m <- glm(Q.kplus1 ~ offset(off), data = data.frame(Q.kplus1, 
                                                             h, off), subset = uncensored & !deterministic, 
                   family = "quasibinomial", weights = h)
          Q1W <- plogis(off + coef(m))
        }
        else {
          m <- glm(Q.kplus1 ~ -1 + h + offset(off), data = data.frame(Q.kplus1, 
                                                                      h, off), subset = uncensored & !deterministic, 
                   family = "quasibinomial")
          Q1W <- plogis(off + coef(m) * h)
        }
      }
      Q1W[deterministic] <- 1
      return(Q1W)
    }
    estg <- function(d, form, Inodes, Ynodes) {
      n <- nrow(d)
      n.g <- length(form)
      gmat <- matrix(NA, nrow = n, ncol = n.g)
      uncensored <- rep(TRUE, n)
      deterministic <- rep(FALSE, n)
      for (i in 1:n.g) {
        if (any(Ynodes < Inodes[i])) {
          Ynode.prev <- max(Ynodes[Ynodes < Inodes[i]])
          deterministic <- d[, Ynode.prev] == 1
        }
        m <- glm(as.formula(form[i]), data = d, subset = uncensored & 
                   !deterministic, family = "binomial")
        gmat[, i] <- predict(m, newdata = d, type = "response")
        gmat[deterministic, i] <- 1
        uncensored <- d[, Inodes[i]] == 1
      }
      return(gmat)
    }
    Bound <- function(x, bounds) {
      x[x < min(bounds)] <- min(bounds)
      x[x > max(bounds)] <- max(bounds)
      return(x)
    }
    n <- nrow(d)
    n.Q <- length(Lnodes)
    n.g <- length(Inodes)
    g1W <- estg(d, gform, Inodes, Ynodes)
    cum.g1W <- Bound(t(apply(g1W, 1, cumprod)), c(gbd, 1))
    empirical.meanwt <- mean(1/cum.g1W[, n.g], na.rm = TRUE)
    cum.g1W[is.na(cum.g1W)] <- Inf
    iptw <- mean(d[, Lnodes[n.Q]] * d[, Inodes[n.g]]/cum.g1W[, 
                                                             n.g])
    wt.n <- 1/cum.g1W[, n.g]/empirical.meanwt
    iptw.wt.n <- mean(d[, Lnodes[n.Q]] * d[, Inodes[n.g]] * wt.n)
    Qstar <- Qinit <- d[, Lnodes[n.Q]]
    IC <- rep(0, n)
    for (i in n.Q:1) {
      Inode.cur <- which.max(Inodes[Inodes < Lnodes[i]])
      uncensored <- d[, Inodes[Inode.cur]] == 1
      if (any(Ynodes < Lnodes[i])) {
        Ynode.prev <- max(Ynodes[Ynodes < Lnodes[i]])
        deterministic <- d[, Ynode.prev] == 1
      }
      else {
        deterministic <- rep(FALSE, n)
      }
      Qinit <- estQ(Qinit, d, Qform[i], uncensored, deterministic, 
                    family = family)
      Qstar.kplus1 <- Qstar
      Qstar <- estQ(Qstar.kplus1, d, Qform[i], uncensored, 
                    deterministic, h = 1/cum.g1W[, Inode.cur], family = family)
      IC[uncensored] <- (IC + (Qstar.kplus1 - Qstar)/cum.g1W[, 
                                                             Inode.cur])[uncensored]
    }
    IC <- IC + Qstar - mean(Qstar)
    return(c(iptw = iptw, iptw.wt.n = iptw.wt.n, Gcomp = mean(Qinit), 
             tmle = mean(Qstar), var.tmle = var(IC)/n))
  }
LtmleFromInputs <- function (inputs)
{
  msm.result <- LtmleMSMFromInputs(inputs)
  if (inputs$verbose){
    message("LtmleFromInputs: preparing results.")
  }
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes %in% 1:2)
  if (num.regimes == 2) {
    class(msm.result) <- "ltmleEffectMeasures"
    return(msm.result)
  }
  names(msm.result$beta.iptw) <- names(msm.result$beta) <- NULL
  iptw <- plogis(msm.result$beta.iptw)
  iptw.list <- list(iptw.estimate = iptw, iptw.IC = iptw *
                      (1 - iptw) * msm.result$IC.iptw[, 1])
  r <- list()
  if (inputs$iptw.only) {
    tmle <- NA
    tmle.IC <- rep(NA, nrow(inputs$data))
  }
  else {
    tmle <- plogis(msm.result$beta)
    tmle.IC <- msm.result$IC[, 1]
  }
  r$estimates <- c(tmle = tmle, iptw = iptw.list$iptw.estimate)
  r$IC <- list(tmle = tmle.IC * tmle * (1 - tmle), iptw = iptw.list$iptw.IC)
  if (!is.null(msm.result$variance.estimate)) {
    stopifnot(length(msm.result$variance.estimate) == 1)
    r$variance.estimate <- msm.result$variance.estimate[1] *
      (tmle * (1 - tmle))^2
  }
  if (inputs$gcomp) {
    names(r$estimates)[1] <- names(r$IC)[1] <- "gcomp"
  }
  r$cum.g <- AsMatrix(msm.result$cum.g[, , 1])
  r$cum.g.unbounded <- AsMatrix(msm.result$cum.g.unbounded[,
                                                           , 1])
  r$cum.g.used <- AsMatrix(msm.result$cum.g.used[, , 1])
  r$gcomp <- inputs$gcomp
  r$fit <- msm.result$fit
  r$fit$g <- r$fit$g[[1]]
  r$fit$Q <- r$fit$Q[[1]]
  r$Qstar <- msm.result$Qstar[, 1, 1]
  r$formulas <- msm.result$formulas
  r$binaryOutcome <- msm.result$binaryOutcome
  r$transformOutcome <- msm.result$transformOutcome == TRUE
  if (msm.result$transformOutcome) {
    Yrange <- attr(msm.result$transformOutcome, "Yrange")
    r$estimates <- r$estimates * diff(Yrange) + min(Yrange)
    r$IC <- lapply(r$IC, function(IC) IC * diff(Yrange))
    r$variance.estimate <- r$variance.estimate * (diff(Yrange))^2
  }
  class(r) <- "ltmle"
  return(r)
}
ltmleMSM <-
  function (data, Anodes, Cnodes = NULL, Lnodes = NULL, Ynodes, 
            survivalOutcome = NULL, Qform = NULL, gform = NULL, gbounds = c(0.01, 
                                                                            1), Yrange = NULL, deterministic.g.function = NULL, SL.library = "glm", 
            SL.cvControl = list(), regimes, working.msm, summary.measures, 
            final.Ynodes = NULL, stratify = FALSE, msm.weights = "empirical", 
            estimate.time = TRUE, gcomp = FALSE, iptw.only = FALSE, deterministic.Q.function = NULL, 
            variance.method = "tmle", observation.weights = NULL, id = NULL) 
  {
    data <- CheckData(data)
    inputs <- CreateInputs(data, Anodes, Cnodes, Lnodes, Ynodes, 
                           survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, 
                           SL.library, SL.cvControl, regimes, working.msm, summary.measures, 
                           final.Ynodes, stratify, msm.weights, estimate.time, gcomp, 
                           iptw.only, deterministic.Q.function, variance.method, 
                           observation.weights, id)
    result <- LtmleMSMFromInputs(inputs)
    result$call <- match.call()
    class(result) <- "ltmleMSM"
    return(result)
  }
LtmleMSMFromInputs <-
  function (inputs)
  {
    if (inputs$estimate.time)
      EstimateTime(inputs)
    if (inputs$verbose){ message("LtmleMSMFromInputs: starting main calculations ...")}
    result <- MainCalcs(inputs)
    result$gcomp <- inputs$gcomp
    result$formulas <- list(Qform = inputs$Qform, gform = inputs$gform)
    result$binaryOutcome <- inputs$binaryOutcome
    result$transformOutcome <- inputs$transformOutcome
    result$survivalOutcome <- inputs$survivalOutcome
    return(result)
  }
MainCalcs <- function (inputs){
  num.final.Ynodes <- length(inputs$final.Ynodes)
  num.betas <- dim(inputs$combined.summary.measures)[2]
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  Qstar <- array(dim = c(n, num.regimes, num.final.Ynodes))
  if (inputs$verbose){ message("MainCalcs: getting weigths ...")}
  all.msm.weights <- GetMsmWeights(inputs)
  new.var.y <- array(dim = c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)
  IC.y <- array(dim = c(n, num.betas, num.final.Ynodes))
  if (inputs$verbose){ message("MainCalcs: estimating G ...")}
  g.list <- EstimateG(inputs)
  if (inputs$verbose){ message("MainCalcs: calculating IPTW ...")}
  iptw <- CalcIPTW(inputs, g.list$cum.g, all.msm.weights)
  fit <- list(g = g.list$fit)
  if (inputs$iptw.only) {
    beta <- rep(NA, length(iptw$beta))
    fitted.msm <- NULL
    variance.estimate <- NULL
    fixed.tmle <- list(cum.g.used = array(NA, dim = dim(g.list$cum.g)))
  }
  else {
    if (inputs$verbose){ message("MainCalcs: calculating TMLE update of Q.")}
    for (j in 1:num.final.Ynodes) {
      fixed.tmle <- FixedTimeTMLE(inputs, nodes = SubsetNodes(inputs$all.nodes,
                                                              final.Ynode = inputs$final.Ynodes[j]), msm.weights = drop3(all.msm.weights[,
                                                                                                                                         , j, drop = FALSE]), combined.summary.measures = dropn(inputs$combined.summary.measures[,
                                                                                                                                                                                                                                 , , j, drop = FALSE], n = 4), g.list = g.list)
      IC <- IC + fixed.tmle$IC
      IC.y[, , j] <- fixed.tmle$IC
      Qstar[, , j] <- fixed.tmle$Qstar
      new.var.y[, , j] <- fixed.tmle$est.var
    }
    fit <- c(fit, fixed.tmle$fit)
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance",
                    exact = TRUE))) {
      return(list(IC = matrix(NA, 1, 1), msm = NULL, beta = qlogis(mean(Qstar)),
                  cum.g = g.list$cum.g, cum.g.unbounded = g.list$cum.g.unbounded,
                  fit = fit, variance.estimate = NULL, beta.iptw = iptw$beta,
                  IC.iptw = iptw$IC, Qstar = Qstar, cum.g.used = fixed.tmle$cum.g.used))
    }
    ## if (inputs$verbose){ message("MainCalcs: fitting pooled MSM.")}
    fitted.msm <- FitPooledMSM(inputs$working.msm, Qstar,
                               inputs$combined.summary.measures, all.msm.weights *
                                 inputs$observation.weights)
    IC <- FinalizeIC(IC, inputs$combined.summary.measures,
                     Qstar, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights,
                     inputs$id)
    C.old <- NormalizeIC(IC, inputs$combined.summary.measures,
                         fitted.msm$m.beta, all.msm.weights, inputs$observation.weights,
                         g.ratio = NULL)
    g.ratio <- CalcGUnboundedToBoundedRatio(g.list, inputs$all.nodes,
                                            inputs$final.Ynodes)
    CheckForVarianceWarning(inputs, g.ratio)
    if (inputs$variance.method == "ic") {
      variance.estimate <- NULL
    }
    else {
      new.var <- matrix(NA, num.betas, num.betas)
      for (i in 1:num.betas) {
        for (j in 1:num.betas) {
          if (num.final.Ynodes > 1) {
            cov.IC <- cov(IC.y[, i, ], IC.y[, j, ])
            diag(cov.IC) <- new.var.y[i, j, ]
            new.var[i, j] <- sum(cov.IC)
          }
          else {
            new.var[i, j] <- new.var.y[i, j, 1]
          }
        }
      }
      C <- NormalizeIC(IC, inputs$combined.summary.measures,
                       fitted.msm$m.beta, all.msm.weights, inputs$observation.weights,
                       g.ratio)
      variance.estimate <- safe.solve(C) %*% new.var %*%
        t(safe.solve(C))
    }
    IC <- t(safe.solve(C.old, t(IC)))
    beta <- coef(fitted.msm$m)
    names(beta) <- inputs$beta.names
  }
  return(list(IC = IC, msm = fitted.msm$m, beta = beta, cum.g = g.list$cum.g,
              cum.g.unbounded = g.list$cum.g.unbounded, fit = fit,
              variance.estimate = variance.estimate, beta.iptw = iptw$beta,
              IC.iptw = iptw$IC, Qstar = Qstar, cum.g.used = fixed.tmle$cum.g.used))
}
MaintainControl <-
  function (data, current.node, nodes) 
  {
    Anodes <- nodes$A
    if (!(current.node %in% Anodes)) 
      return(NULL)
    if (!(any(Anodes < current.node))) 
      return(NULL)
    prev.a.node <- max(Anodes[Anodes < current.node])
    is.deterministic <- data[, prev.a.node] == 0
    return(list(is.deterministic = is.deterministic, prob1 = 0))
  }
MaintainTreatment <-
  function (data, current.node, nodes) 
  {
    Anodes <- nodes$A
    if (!(current.node %in% Anodes)) 
      return(NULL)
    if (!(any(Anodes < current.node))) 
      return(NULL)
    prev.a.node <- max(Anodes[Anodes < current.node])
    is.deterministic <- data[, prev.a.node] == 1
    return(list(is.deterministic = is.deterministic, prob1 = 1))
  }
MakePosDef <-
  function (variance.estimate) 
  {
    orig.variance.estimate <- variance.estimate
    try.result <- try({
      near.pd <- Matrix::nearPD(variance.estimate)
      variance.estimate <- as.matrix(near.pd$mat)
    }, silent = TRUE)
    if (inherits(try.result, "try-error") || !near.pd$converged || 
        any(abs(orig.variance.estimate - variance.estimate) > 
            0.001 & (abs(orig.variance.estimate - variance.estimate)/orig.variance.estimate) > 
            0.10000000000000001)) {
      warning("Covariance matrix from EstimateVariance not positive definite, unable to compute standard errors. You may want to try variance.method='ic'.")
      variance.estimate <- matrix(nrow = nrow(variance.estimate), 
                                  ncol = ncol(variance.estimate))
    }
    return(variance.estimate)
  }
Mode <-
  function (x, na.rm = FALSE) 
  {
    if (na.rm) 
      x <- x[!is.na(x)]
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
NodeToIndex <-
  function (data, node) 
  {
    if (is.numeric(node) || is.null(node)) 
      return(node)
    if (!is.character(node)) 
      stop("nodes must be numeric, character, or NULL")
    index <- match(node, names(data))
    if (anyNA(index)) {
      stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
    }
    return(index)
  }
NormalizeIC <-
  function (IC, combined.summary.measures, m.beta, msm.weights, 
            observation.weights, g.ratio) 
  {
    n <- dim(combined.summary.measures)[1]
    num.betas <- dim(combined.summary.measures)[2]
    num.regimes <- dim(combined.summary.measures)[3]
    num.final.Ynodes <- dim(combined.summary.measures)[4]
    if (is.null(g.ratio)) {
      g.ratio <- array(1, dim = c(n, num.regimes, num.final.Ynodes))
    }
    C <- array(0, dim = c(num.betas, num.betas))
    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimes) {
        tempC <- crossprod(combined.summary.measures[, , 
                                                     i, j] * g.ratio[, i, j], combined.summary.measures[, 
                                                                                                        , i, j] * g.ratio[, i, j] * msm.weights[, i, 
                                                                                                                                                j] * m.beta[, i, j] * (1 - m.beta[, i, j]) * 
                             observation.weights)
        if (anyNA(tempC)) 
          stop("NA in tempC")
        C <- C + tempC
      }
    }
    C <- C/n
    if (rcond(C) < 9.9999999999999998e-13) {
      C <- matrix(NA, nrow = num.betas, ncol = num.betas)
      warning("rcond(C) near 0, standard errors not available")
    }
    return(C)
  }
print.ltmle <-
  function (x, ...) 
  {
    PrintCall(x$call)
    if (x$gcomp) {
      cat("GCOMP Estimate: ", x$estimates["gcomp"], "\n")
    }
    else {
      cat("TMLE Estimate: ", x$estimates["tmle"], "\n")
    }
    invisible(x)
  }
print.ltmleEffectMeasures <-
  function (x, ...) 
  {
    PrintCall(x$call)
    cat("Use summary(...) to get estimates, standard errors, p-values, and confidence intervals for treatment EYd, control EYd, additive effect, relative risk, and odds ratio.\n")
    invisible(x)
  }
print.ltmleMSM <-
  function (x, ...) 
  {
    PrintCall(x$call)
    if (x$gcomp) {
      cat("GCOMP Beta Estimates: \n")
    }
    else {
      cat("TMLE Beta Estimates: \n")
    }
    print(x$beta)
    if (x$transformOutcome) {
      Yrange <- attr(x$transformOutcome, "Yrange")
      cat("NOTE: The MSM is modeling the transformed outcome ( Y -", 
          min(Yrange), ")/(", max(Yrange), "-", min(Yrange), 
          ")")
    }
    invisible(x)
  }
print.summary.ltmle <-
  function (x, ...) 
  {
    cat("Estimator: ", x$estimator, "\n")
    if (x$estimator == "gcomp") {
      cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")
    }
    PrintCall(x$call)
    PrintSummary(x$treatment)
    CheckVarianceEstimateRatio(x)
    invisible(x)
  }
print.summary.ltmleEffectMeasures <-
  function (x, ...) 
  {
    cat("Estimator: ", x$estimator, "\n")
    if (x$estimator == "gcomp") {
      cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")
    }
    PrintCall(x$call)
    lapply(x$effect.measures, PrintSummary)
    CheckVarianceEstimateRatio(x)
    invisible(x)
  }
print.summary.ltmleMSM <-
  function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), 
            ...) 
  {
    cat("Estimator: ", x$estimator, "\n")
    if (x$estimator == "gcomp") {
      cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")
    }
    printCoefmat(x$cmat, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", has.Pvalue = TRUE, ...)
    if (x$transformOutcome) {
      Yrange <- attr(x$transformOutcome, "Yrange")
      cat("NOTE: The MSM is modeling the transformed outcome ( Y -", 
          min(Yrange), ")/(", max(Yrange), "-", min(Yrange), 
          ")")
    }
    CheckVarianceEstimateRatio(x)
    invisible(x)
  }
PrintCall <-
  function (cl) 
  {
    cat("Call:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  }
PrintSummary <-
  function (x) 
  {
    if (!is.null(x$long.name)) 
      cat(x$long.name, ":\n", sep = "")
    cat("   Parameter Estimate: ", signif(x$estimate, 5), "\n")
    if (x$log.std.err) {
      if (x$long.name == "Relative Risk") {
        param.abbrev <- "RR"
      }
      else if (x$long.name == "Odds Ratio") {
        param.abbrev <- "OR"
      }
      else {
        stop("unexpected x$long.name")
      }
      cat("  Est Std Err log(", param.abbrev, "):  ", sep = "")
    }
    else {
      cat("    Estimated Std Err:  ")
    }
    cat(signif(x$std.dev, 5), "\n")
    cat("              p-value: ", ifelse(x$pvalue <= 2 * 10^-16, 
                                          "<2e-16", signif(x$pvalue, 5)), "\n")
    cat("    95% Conf Interval:", paste("(", signif(x$CI[1], 
                                                    5), ", ", signif(x$CI[2], 5), ")", sep = ""), "\n\n")
    invisible(x)
  }
RegimesFromAbar <-
  function (data, abar, rule) 
  {
    if (!is.null(rule)) {
      if (!(missing(abar) || is.null(abar))) 
        stop("'abar' should not be specified when using a 'rule' function")
      rule.output.length <- length(as.numeric(rule(data[1, 
      ])))
      if (all(sapply(data, is.numeric))) {
        abar <- apply(data, 1, rule)
        if (rule.output.length == 1) {
          abar <- AsMatrix(abar)
        }
        else {
          abar <- t(abar)
        }
      }
      else {
        if (nrow(data) > 10000) 
          warning("Using a rule input may be significantly slower than using abar/regimes if your data contains censoring nodes or is otherwise not all numeric.")
        abar <- matrix(nrow = nrow(data), ncol = rule.output.length)
        for (i in 1:nrow(data)) {
          abar[i, ] <- as.numeric(rule(data[i, ]))
        }
      }
    }
    if (is.vector(abar)) {
      abar <- matrix(rep(abar, each = nrow(data)), nrow = nrow(data))
    }
    else if (is.null(abar)) {
      abar <- matrix(numeric(0), nrow = nrow(data), ncol = 0)
    }
    else if (is.function(abar)) {
      stop("abar should be a vector or matrix, not a function. Use the 'rule' paramter to pass a function instead.")
    }
    regimes <- abar
    dim(regimes) <- c(nrow(regimes), ncol(regimes), 1)
    return(regimes)
  }
ReorderFits <-
  function (l1) 
  {
    if (length(l1) == 0) 
      l1 <- list("no fit due to no A/C nodes")
    num.regimes <- length(l1[[1]])
    num.nodes <- length(l1)
    l2 <- vector("list", num.regimes)
    for (i in 1:num.regimes) {
      l2[[i]] <- vector("list", num.nodes)
      names(l2[[i]]) <- names(l1)
      for (j in 1:num.nodes) {
        l2[[i]][[j]] <- l1[[j]][[i]]
      }
    }
    return(l2)
  }
repmat <-
  function (X, m, n) 
  {
    mx <- dim(X)[1]
    nx <- dim(X)[2]
    if ((m == 0) || (n == 0)) 
      return(matrix(numeric(0), nrow = mx * m, ncol = nx * 
                      n))
    return(matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T))
  }
rexpit <-
  function (x) 
    rbinom(n = length(x), size = 1, prob = plogis(x))
RhsVars <-
  function (f) 
  {
    f <- as.formula(f)
    return(all.vars(f[[3]]))
  }
safe.solve <-
  function (a, b) 
  {
    if (missing(b)) {
      try.result <- try(x <- solve(a))
    }
    else {
      try.result <- try(x <- solve(a, b))
    }
    if (inherits(try.result, "try-error")) {
      if (missing(b)) {
        x <- matrix(nrow = nrow(a), ncol = ncol(a))
      }
      else {
        x <- matrix(nrow = ncol(a), ncol = ncol(AsMatrix(b)))
      }
      warning("Error in solve(), standard errors not available")
    }
    return(x)
  }
Scale <-
  function (x, min.y, max.y) 
  {
    if (all(is.na(x))) 
      stop("all NA in Scale")
    r <- range(x, na.rm = TRUE)
    if (diff(r) > 0) {
      return((x - r[1])/diff(r) * (max.y - min.y) + min.y)
    }
    else {
      if (r[1] >= min.y && r[1] <= max.y) {
        return(rep(r[1], length(x)))
      }
      else {
        return(rep(mean(c(min.y, max.y)), length(x)))
      }
    }
  }
SetA <-
  function (data, regimes, Anodes, regime.index, cur.node) 
  {
    Anode.index <- which(Anodes < cur.node)
    for (a in Anode.index) {
      data[, Anodes[a]] <- regimes[, a, regime.index]
    }
    return(data)
  }
SetSeedIfRegressionTesting <-
  function () 
  {
    seed <- Sys.getenv("LTMLE.REGRESSION.TESTING.SEED")
    stopifnot(length(seed) == 1)
    if (seed != "") {
      seed <- as.numeric(seed)
      stopifnot(is.finite(seed))
      set.seed(seed)
    }
    invisible(NULL)
  }
ShowGlmMessage <-
  function () 
  {
    message("speedglm failed, using glm instead. If you see a lot of this message and you have large absolute values in data[, Lnodes], you may get a speed performance improvement by rescaling these values.")
  }
sseq <-
  function (from, to) 
  {
    if (from > to) 
      return(integer(0))
    seq(from, to)
  }
sub2ind <-
  function (row, col, num.rows) 
  {
    return((col - 1) * num.rows + row)
  }
SubsetNodes <-
  function (nodes, final.Ynode) 
  {
    return(lapply(nodes, function(x) x[x <= final.Ynode]))
  }
summarize_analysis <- function(x,outcome,time=1:5,outcome_list,raw=FALSE,merge=FALSE){
  x=copy(x)
  if (length(time)>0){
    if (match(outcome,names(x),nomatch=0))
      x=rbindlist(x[[outcome]][time])
    else # for subgroup analyses
      x=rbindlist(x[time])
  } else{
    x=x[[outcome]]
  }
  x[what!="rr",estimate:=100*estimate]
  x[what!="rr",lower:=100*lower]
  x[what!="rr",upper:=100*upper]
  x[,treatment:=factor(treatment,
                       levels=c("glp1_and_other","glp1_and_sglt2","other_and_other","sglt2_inhib_and_other"),
                       labels=c("GLP1-RA and DPP4/SU/TZD","GLP1-RA and SGLT2i","Dual DPP4/SU/TZD","SGLT2i and DPP4/SU/TZD"))]
  ## x[,reference:=factor(reference,levels=c("","glp1_and_sglt2"),labels=c("","GLP1-RA and SGLT2i"))]
  x[,reference:=factor(reference,levels=c("","other_and_other"),labels=c("","Dual DPP4/SU/TZD"))]
  if (raw){
    if ("subset"%in%names(x)){
      R=x[reference==""][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"estimate"=estimate,lower=lower,upper=upper)]
      D=x[what=="ate"&reference!=""][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"estimate"=estimate,lower=lower,upper=upper,"p-value"=format.pval(2*pnorm(abs(estimate/100)/se,lower.tail=FALSE),eps=0.001,digits=2))]
      RR=x[what=="rr"&reference!=""][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"estimate"=estimate,lower=lower,upper=upper,"p-value"="")]
    }else{
      R=x[reference==""][,.("Years"=horizon,"Regimen"=treatment,N=N,"estimate"=estimate,lower=lower,upper=upper)]
      D=x[what=="ate"&reference!=""][,.("Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"estimate"=estimate,lower=lower,upper=upper,"p-value"=format.pval(2*pnorm(abs(estimate/100)/se,lower.tail=FALSE),eps=0.001,digits=2))]
      RR=x[what=="rr"&reference!=""][,.("Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"estimate"=estimate,lower=lower,upper=upper,"p-value"="")]
    }
  }else{
    if ("subset"%in%names(x)){
      R=x[reference==""][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"Risk"=formatCI(x=estimate,lower=lower,upper=upper,digits=1,show.x=TRUE))]
      D=x[what=="ate"&(is.na(reference)|reference!="")][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"Difference"=formatCI(x=estimate,lower=lower,upper=upper,show.x=TRUE,digits=2),"p-value"=format.pval(2*pnorm(abs(estimate/100)/se,lower.tail=FALSE),eps=0.001,digits=2))]
      RR=x[what=="rr"&(is.na(reference)|reference!="")][,.("Subset"=subset,"Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"Ratio"=formatCI(x=estimate,lower=lower,upper=upper,show.x=TRUE,digits=2),"p-value"="")]
    }else{
      R=x[reference==""][,.("Years"=horizon,"Regimen"=treatment,N=N,"Risk"=formatCI(x=estimate,lower=lower,upper=upper,digits=1,show.x=TRUE))]
      D=x[what=="ate"&(is.na(reference)|reference!="")][,.("Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"Difference"=formatCI(x=estimate,lower=lower,upper=upper,show.x=TRUE,digits=2),
                                                           "p-value"=format.pval(2*pnorm(abs(estimate/100)/se,lower.tail=FALSE),eps=0.001,digits=2))]
      RR=x[what=="rr"&(is.na(reference)|reference!="")][,.("Years"=horizon,"Regimen"=treatment,N=N,"Reference"=reference,"Ratio"=formatCI(x=estimate,lower=lower,upper=upper,show.x=TRUE,digits=2),"p-value"="")]
      
    }
  }
  if (match(outcome,names(outcome_list),nomatch=0)){
    out <- list(R=cbind(Outcome=outcome_list[[outcome]],R),
                D=cbind(Outcome=outcome_list[[outcome]],D),
                RR=cbind(Outcome=outcome_list[[outcome]],RR))
  } else{
    out <- list(R=R,D=D,RR=RR)
  }
  if (merge){
    out <- merge(R,D,all.y=TRUE,all.x=TRUE,on=intersect(names(R),names(D)))
    out[,Regimen:=factor(Regimen,levels=c("GLP1-RA and SGLT2i","GLP1-RA and DPP4/SU/TZD","SGLT2i and DPP4/SU/TZD","Dual DPP4/SU/TZD"))]
    y <- out[["Years"]][1]
    out[,Years:=NULL]
    out[is.na(Reference),Reference:=""]
    out[is.na(Difference)][["p-value"]] <- ""
    out[is.na(Difference),Difference:="Reference"]
    out[,Reference:=NULL]
    setnames(out,"Risk",paste0(y,"-year  (%)risk"))
    setnames(out,"Difference",paste0(y,"-year risk difference (%)"))
    if ("Subset"%in%names(out)){
      setkey(out,Subset,Regimen)
      setkey(RR,Subset,Regimen)
      corder <- names(out)
      out <- RR[,.(Subset,Regimen,"Ratio"=Ratio)][out]
    }else{
      setkey(out,Regimen)
      setkey(RR,Regimen)
      corder <- names(out)
      out <- RR[,.(Regimen,"Ratio"=Ratio)][out]
    }
    setcolorder(out,c(corder,"Ratio"))
    out[Regimen=="Dual DPP4/SU/TZD",Ratio:="Reference"]
    setnames(out,"Ratio",paste0(y,"-year risk ratio"))
  }
  out[]
}
summary.ltmle <-
  function (object, estimator = ifelse(object$gcomp, "gcomp", "tmle"),
            ...)
  {
    if ("control.object" %in% names(list(...)))
      stop("control.object has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
    if (!estimator[1] %in% c("tmle", "iptw", "gcomp"))
      stop("estimator should be one of: tmle, iptw, gcomp. If you are trying to use control.object, the control.object parameter has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
    if (estimator == "tmle" && object$gcomp)
      stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
    if (estimator == "gcomp" && !object$gcomp)
      stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
    IC.variance <- var(object$IC[[estimator]])
    if (estimator == "tmle" && !is.null(object$variance.estimate)) {
      v <- max(IC.variance, object$variance.estimate)
    }
    else {
      v <- IC.variance
    }
    variance.estimate.ratio = v/IC.variance
    if (object$binaryOutcome) {
      CIBounds <- c(0, 1)
    }
    else {
      CIBounds <- c(-Inf, Inf)
    }
    treatment <- GetSummary(list(long.name = NULL, est = object$estimates[estimator],
                                 gradient = 1, log.std.err = FALSE, CIBounds = CIBounds),
                            v,
                            n = length(object$IC[[estimator]]))
    ans <- list(treatment = treatment, call = object$call, estimator = estimator,
                variance.estimate.ratio = variance.estimate.ratio)
    class(ans) <- "summary.ltmle"
    return(ans)
  }
summary.ltmleEffectMeasures <-
  function (object, estimator = ifelse(object$gcomp, "gcomp", "tmle"), 
            ...) 
  {
    info <- GetSummaryLtmleMSMInfo(object, estimator)
    beta <- info$estimate
    IC <- info$IC
    y0 <- plogis(beta[1])
    y1 <- plogis(beta[1] + beta[2])
    names(y0) <- names(y1) <- NULL
    eff.list <- list(treatment = list(long.name = "Treatment Estimate", 
                                      est = y1, gradient = c(y1 * (1 - y1), y1 * (1 - y1)), 
                                      log.std.err = FALSE, CIBounds = 0:1), control = list(long.name = "Control Estimate", 
                                                                                           est = y0, gradient = c(y0 * (1 - y0), 0), log.std.err = FALSE, 
                                                                                           CIBounds = 0:1), ATE = list(long.name = "Additive Treatment Effect", 
                                                                                                                       est = y1 - y0, gradient = c(y1 * (1 - y1) - y0 * (1 - 
                                                                                                                                                                           y0), y1 * (1 - y1)), log.std.err = FALSE, CIBounds = c(-1, 
                                                                                                                                                                                                                                  1)), RR = list(long.name = "Relative Risk", est = y1/y0, 
                                                                                                                                                                                                                                                 gradient = c(y0 - y1, 1 - y1), log.std.err = TRUE, CIBounds = c(0, 
                                                                                                                                                                                                                                                                                                                 Inf)), OR = list(long.name = "Odds Ratio", est = exp(beta[2]), 
                                                                                                                                                                                                                                                                                                                                  gradient = c(0, 1), log.std.err = TRUE, CIBounds = c(0, 
                                                                                                                                                                                                                                                                                                                                                                                       Inf)))
    if (!object$binaryOutcome) {
      eff.list$RR <- eff.list$OR <- NULL
    }
    n <- nrow(IC)
    measures.IC <- lapply(eff.list, GetSummary, var(IC), n)
    if (is.null(object$variance.estimate)) {
      measures.variance.estimate <- NULL
    }
    else {
      measures.variance.estimate <- lapply(eff.list, GetSummary, 
                                           object$variance.estimate, n)
    }
    measures.max <- measures.IC
    for (i in seq_along(measures.variance.estimate)) {
      std.dev.diff <- measures.variance.estimate[[i]]$std.dev - 
        measures.IC[[i]]$std.dev
      if (!is.na(std.dev.diff) && (std.dev.diff > 0)) {
        measures.max[[i]] <- measures.variance.estimate[[i]]
      }
    }
    if (object$transformOutcome) {
      Yrange <- attr(object$transformOutcome, "Yrange")
      measures.max <- lapply(measures.max, function(x) {
        x$estimate <- x$estimate * diff(Yrange)
        x$std.dev <- x$std.dev * diff(Yrange)
        x$CI <- x$CI * diff(Yrange)
        if (x$long.name %in% c("Treatment Estimate", "Control Estimate")) {
          x$estimate <- x$estimate + min(Yrange)
          x$CI <- x$CI + min(Yrange)
        }
        else {
          stopifnot(x$long.name == "Additive Treatment Effect")
        }
        return(x)
      })
    }
    ans <- list(call = object$call, effect.measures = measures.max, 
                variance.estimate.ratio = info$variance.estimate.ratio, 
                estimator = estimator)
    class(ans) <- "summary.ltmleEffectMeasures"
    return(ans)
  }
summary.ltmleMSM <-
  function (object, estimator = ifelse(object$gcomp, "gcomp", "tmle"), 
            ...) 
  {
    info <- GetSummaryLtmleMSMInfo(object, estimator)
    estimate <- info$estimate
    v <- info$v
    n <- nrow(info$IC)
    std.dev <- sqrt(v/n)
    pval <- 2 * GetPValue(-abs(estimate/std.dev), n)
    CI <- GetCI(estimate, std.dev, n)
    cmat <- cbind(estimate, std.dev, CI, pval)
    dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", 
                                              "CI 2.5%", "CI 97.5%", "p-value"))
    ans <- list(cmat = cmat, estimator = estimator, transformOutcome = object$transformOutcome, 
                variance.estimate.ratio = info$variance.estimate.ratio)
    class(ans) <- "summary.ltmleMSM"
    return(ans)
  }
SuppressGivenWarnings <-
  function (expr, warningsToIgnore) 
  {
    h <- function(w) {
      if (w$message %in% warningsToIgnore) 
        invokeRestart("muffleWarning")
    }
    withCallingHandlers(expr, warning = h)
  }
TransformOutcomes <-
  function (data, nodes, Yrange) 
  {
    all.Y <- unlist(data[, nodes$Y])
    transformOutcome <- FALSE
    if (!is.null(Yrange)) {
      rng <- range(all.Y, na.rm = TRUE)
      if (min(rng) < min(Yrange) || max(rng) > max(Yrange)) {
        message("Some Ynodes are not in [Yrange[1], Yrange[2]], Y values are truncated")
        data[, nodes$Y][data[, nodes$Y] < min(Yrange)] <- min(Yrange)
        data[, nodes$Y][data[, nodes$Y] > max(Yrange)] <- max(Yrange)
      }
      transformOutcome <- TRUE
    }
    else {
      Yrange <- range(all.Y, na.rm = TRUE)
      if (min(Yrange) < 0 || max(Yrange) > 1) {
        transformOutcome <- TRUE
        message("Some Ynodes are not in [0, 1], and Yrange was NULL, so all Y nodes are being\ntransformed to (Y-min.of.all.Ys)/range.of.all.Ys")
      }
    }
    if (transformOutcome) {
      attr(transformOutcome, "Yrange") <- Yrange
      data[, nodes$Y] <- (data[, nodes$Y] - min(Yrange))/diff(Yrange)
    }
    return(list(data = data, transformOutcome = transformOutcome))
  }
UpdateQ <-
  function (Qstar.kplus1, logitQ, combined.summary.measures, cum.g, 
            working.msm, uncensored, intervention.match, is.deterministic, 
            msm.weights, gcomp, observation.weights) 
  {
    n <- nrow(logitQ)
    num.regimes <- ncol(logitQ)
    off <- as.vector(logitQ)
    Y <- as.vector(Qstar.kplus1)
    if (length(cum.g) == 0) 
      cum.g <- 1
    stacked.summary.measures <- apply(combined.summary.measures, 
                                      2, rbind)
    subs.vec <- uncensored & !is.deterministic & as.vector(intervention.match)
    ## print(head(cum.g))
    weight.vec <- numeric(n * num.regimes)
    weight.vec[subs.vec] <- (observation.weights * as.vector(msm.weights)/as.vector(cum.g))[subs.vec]
    if (anyNA(weight.vec)) 
      stop("NA in weight.vec")
    f <- as.formula(paste(working.msm, "+ offset(off)"))
    data.temp <- data.frame(Y, stacked.summary.measures, off)
    if (gcomp) {
      Qstar <- plogis(logitQ)
      m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
    }
    else {
      if (any(weight.vec > 0)) {
        m <- ltmle.glm(f, data = data.temp[weight.vec > 0, 
        ], family = quasibinomial(), weights = as.vector(scale(weight.vec[weight.vec > 
                                                                            0], center = FALSE)))
        Qstar <- matrix(predict(m, newdata = data.temp, type = "response"), 
                        nrow = nrow(logitQ))
      }
      else {
        Qstar <- plogis(logitQ)
        m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
      }
    }
    indicator <- matrix(uncensored * observation.weights, nrow = nrow(stacked.summary.measures), 
                        ncol = ncol(stacked.summary.measures)) * matrix(intervention.match, 
                                                                        nrow = nrow(stacked.summary.measures), ncol = ncol(stacked.summary.measures))
    h.g.ratio <- stacked.summary.measures/matrix(cum.g, nrow = nrow(stacked.summary.measures), 
                                                 ncol = ncol(stacked.summary.measures)) * indicator
    dim(h.g.ratio) <- c(n, num.regimes, ncol(h.g.ratio))
    for (i in 1:num.regimes) {
      h.g.ratio[, i, ] <- h.g.ratio[, i, ] * msm.weights[, 
                                                         i]
      weight.zero.index <- msm.weights[, i] == 0
      h.g.ratio[weight.zero.index, i, ] <- 0
    }
    cum.g.used <- weight.vec > 0 & msm.weights > 0
    dim(cum.g.used) <- c(n, num.regimes)
    return(list(Qstar = Qstar, h.g.ratio = h.g.ratio, X = stacked.summary.measures, 
                off = off, fit = m, cum.g.used = cum.g.used))
  }
VarianceAvailableWarning <-
  function (inputs) 
  {
    if (!inputs$binaryOutcome) 
      return("Robust variance estimate is not currently available with non binary outcomes")
    if (!is.null(inputs$deterministic.Q.function)) 
      return("Robust variance estimate is not currently available with deterministic.Q.function")
    if (inputs$gcomp) 
      return("Robust variance estimate is not currently available with gcomp")
    if (inputs$stratify) 
      return("Robust variance estimate is not currently available with stratify=TRUE")
    if (!is.null(inputs$id)) 
      return("Robust variance estimate is not currently available with non-null id")
    return(NULL)
  }


event_node_manipulator <- function(data,k,outcome,competing="Dead",censored="Censored",outcome_is_competing="death"){
  #
  # manipulation of the event nodes
  #
  Y_nodes <- paste0(outcome,"_",1:k)
  D_nodes <- paste0(competing,"_",1:k)
  C_nodes <- paste0(censored,"_",1:k)
  Y_nodes_position <- match(Y_nodes,names(data))
  C_nodes_position <- match(C_nodes,names(data))
  D_nodes_position <- match(D_nodes,names(data))
  #
  # IMPORTANT: when outcome or death occurs in an interval followed by censoring,
  # we choose to uncensor the outcome and the death.
  #
  for (q in 1:k){
    if (q<k){
      has_outcome_or_death_and_censored <- (((data[[Y_nodes_position[[q]]]]%in%1)|(data[[D_nodes_position[[q]]]]%in%1)) &(data[[C_nodes_position[[q]]]]%in%"censored"))
    } else{
      has_outcome_or_death_and_censored <- ((data[[Y_nodes_position[[q]]]]%in%1)&(data[[C_nodes_position[[q]]]]%in%"censored"))
    }
    if (any(has_outcome_or_death_and_censored)){
      set(data,j=C_nodes_position[[q]],i=which(has_outcome_or_death_and_censored),value="uncensored")
    }
  }
  #
  # all nodes (but not outcome and competing risk nodes) should be NA after outcome
  # or competing risk (i.e., death) has occurred
  #
  for (Ok in Y_nodes[-k]){
    later_nodes=setdiff((match(Ok,names(data))+1):NCOL(data),Y_nodes_position)
    if (any(has_outcome <- (data[[Ok]]%in%1))){
      for (l in later_nodes) set(data,j=l,i=which(has_outcome),value=NA)
    }
  }
  if (outcome!=outcome_is_competing){
    # the last Death node occurs after outcome in the last interval and has been removed
    for (Dk in D_nodes[-k]){
      # for some reason Y_nodes must not be missing after death!
      later_nodes=setdiff((match(Dk,names(data))+1):NCOL(data),Y_nodes_position)
      if (any(has_died <- (data[[Dk]]%in%1))){
        for (l in later_nodes) set(data,j=l,i=which(has_died),value=NA)
      }
    }
  }
  #
  # all nodes including outcome should be NA as soon as censored occurred
  #
  for (Ck in C_nodes){
    ## later_nodes=setdiff((match(Ck,names(data))+1):NCOL(data),c(Y_nodes,D_nodes))
    later_nodes=(match(Ck,names(data))+1):NCOL(data)
    if (any(has_censored <- (data[[Ck]]%in%"censored"))){
      for (l in later_nodes) set(data,j=l,i=which(has_censored),value=NA)
    }
  }
  data
}

## Get g-formulas and Q-formulas

get_formulas <- function(time_horizon,
                         work_data,
                         name_outcome,
                         name_baseline_covariates,
                         name_time_covariates,
                         name_regimen,
                         name_censoring = NULL,
                         name_comp.event = NULL,
                         Markov = NULL, ## Names of time varying covariates with Markov property. Note that regimen is assumed NOT to be Markov
                         constant_variables = NULL){
  
  if (length(Markov)>0 && Markov[[1]]!="")
    if (any(not_found <- !(Markov%in%name_time_covariates)))
      stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                  paste(Markov[not_found],collapse=", ")))
  
  time_grid = 0:time_horizon
  K = length(time_grid)
  gform = if(length(name_regimen)==2){
    c(paste0(name_regimen[[1]],"_0"," ~ ", get_rhs(timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                   name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                   regimen = FALSE, Markov = Markov, constant_variables = constant_variables)),
      paste0(name_regimen[[2]],"_0"," ~ ", get_rhs(timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                   name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                   regimen = FALSE, Markov = Markov, constant_variables = constant_variables)),
      if(length(name_censoring)>0){paste0(name_censoring,"_1"," ~ ", get_rhs(timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                                             name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                                             regimen = TRUE, Markov = Markov,
                                                                             constant_variables = constant_variables))} else{})
  } else{
    c(paste0(name_regimen,"_0"," ~ ", get_rhs(timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                              name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                              regimen = FALSE, Markov = Markov, constant_variables = constant_variables)),
      if(length(name_censoring)>0){paste0(name_censoring,"_1"," ~ ", get_rhs(timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                                             name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                                             regimen = TRUE, Markov = Markov, constant_variables = constant_variables))} else{})
  }
  if(time_horizon>1){
    gform = c(gform, unlist(lapply(1:(time_horizon-1),function(tk){
      if(length(name_regimen)==2){
        c(paste0(name_regimen[[1]],"_",tk," ~ ", get_rhs(timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                         name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                         Markov = Markov, constant_variables = constant_variables)),
          paste0(name_regimen[[2]],"_",tk," ~ ", get_rhs(timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                         name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                         Markov = Markov, constant_variables = constant_variables)),
          if(length(name_censoring)>0) {paste0(name_censoring,"_", tk+1, " ~ ",
                                               get_rhs(timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                       name_time_covariates = name_time_covariates,
                                                       name_regimen = name_regimen, regimen = TRUE,
                                                       Markov = Markov, constant_variables = constant_variables))}
          else {})
      } else{
        c(paste0(name_regimen,"_",tk," ~ ", get_rhs(timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                    name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                    Markov = Markov, constant_variables = constant_variables)),
          if(length(name_censoring)>0) {paste0(name_censoring,"_", tk+1, " ~ ",
                                               get_rhs(timepoint = tk+1, name_baseline_covariates = name_baseline_covariates,
                                                       name_time_covariates = name_time_covariates,
                                                       name_regimen = name_regimen, regimen = TRUE,
                                                       Markov = Markov, constant_variables = constant_variables))}
          else {})
      }
    })))
  }
  
  
  ## Note that A_k ~ V + L_0 + ... + L_(k-1) + A_(k-1) for k = 1,..., time_horizon, but A_0 ~ V + L_0
  ## i.e., regimen at baseline depends on additional baseline covariates, whereas in general, regimen depends
  ## on the previously observed covariates and regimen.
  ## The reason for this is we do not want to mistakenly assume that L_1 -> A_1 when in reality A_1 happens before L_1
  
  Qform <- unlist(lapply(1:time_horizon,function(tk){
    paste0("Q.kplus1 ~ ", get_rhs(timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                  name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                  Markov = Markov, constant_variables = constant_variables))
  }))
  
  names(Qform)=paste0(name_outcome,"_",1:time_horizon)
  # names for treatment and censoring formulas
  names(gform) <- as.character(unlist(lapply(gform,function(x){strsplit(x," ~ ")[[1]][[1]]})))
  list(gform = gform, Qform = Qform)
}  

## Get g-formulas and Q-formulas

get_formulas_sim <- function(time_horizon,
                             work_data,
                             name_outcome,
                             name_baseline_covariates,
                             name_time_covariates,
                             name_regimen,
                             name_censoring = NULL,
                             name_comp.event = NULL,
                             Markov = NULL, ## Names of time varying covariates with Markov property. Note that regimen is assumed NOT to be Markov
                             concurrentY, #Whether Yt is predicted from Lt (TRUE, used for sim data from coefficients) or Lt-1 (FALSE)
                             firstL=TRUE, # If L->Y so first block needs to be names by 
                             constant_variables = NULL){
  
  if (length(Markov)>0 && Markov[[1]]!="")
    if (any(not_found <- !(Markov%in%name_time_covariates)))
      stop(paste0("The following variables in argument Markov do not match time_covariates:\n",
                  paste(Markov[not_found],collapse=", ")))
  
  offset=ifelse(concurrentY,0,1) #If not concurrentY, only predict from previous time point
  
  time_grid = 0:time_horizon
  K = length(time_grid)
  gform = if(length(name_regimen)==2){
    c(paste0(name_regimen[[1]],"_0"," ~ ", get_rhs(concurrentY=concurrentY, timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                   name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                   regimen = FALSE, Markov = Markov, constant_variables = constant_variables)),
      paste0(name_regimen[[2]],"_0"," ~ ", get_rhs(concurrentY=concurrentY, timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                                   name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                                   regimen = FALSE, Markov = Markov, constant_variables = constant_variables))#,
      # if(length(name_censoring)>0){paste0(name_censoring,"_1"," ~ ", get_rhs(concurrentY=FALSE, timepoint = 0, name_baseline_covariates = name_baseline_covariates,
      #                                                                        name_time_covariates = name_time_covariates, name_regimen = name_regimen,
      #                                                                        regimen = TRUE, Markov = Markov,
      #                                                                        constant_variables = constant_variables))} else{}
    )
  } else{
    c(paste0(name_regimen,"_0"," ~ ", get_rhs(concurrentY=concurrentY, timepoint = 0, name_baseline_covariates = name_baseline_covariates,
                                              name_time_covariates = name_time_covariates, name_regimen = name_regimen,
                                              regimen = FALSE, Markov = Markov, constant_variables = constant_variables))#,
      # if(length(name_censoring)>0){paste0(name_censoring,"_1"," ~ ", get_rhs(concurrentY=concurrentY, timepoint = 0, name_baseline_covariates = name_baseline_covariates,
      #                                                                        name_time_covariates = name_time_covariates, name_regimen = name_regimen,
      #                                                                        regimen = TRUE, Markov = Markov, constant_variables = constant_variables))} else{}
    )
  }
  if(time_horizon>1){
    gform = c(gform, unlist(lapply(1:(time_horizon-offset),function(tk){
      if(length(name_regimen)==2){
        c(paste0(name_regimen[[1]],"_",tk," ~ ", get_rhs(concurrentY=concurrentY, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                         name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                         Markov = Markov, constant_variables = constant_variables)),
          paste0(name_regimen[[2]],"_",tk," ~ ", get_rhs(concurrentY=concurrentY, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                         name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                         Markov = Markov, constant_variables = constant_variables)),
          if(length(name_censoring)>0) {paste0(name_censoring,"_", tk+offset, " ~ ",
                                               get_rhs(concurrentY=concurrentY, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                       name_time_covariates = name_time_covariates,
                                                       name_regimen = name_regimen, regimen = TRUE,
                                                       Markov = Markov, constant_variables = constant_variables))}
          else {})
      } else{
        c(paste0(name_regimen,"_",tk," ~ ", get_rhs(concurrentY=concurrentY, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                    name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                    Markov = Markov, constant_variables = constant_variables)),
          if(length(name_censoring)>0) {paste0(name_censoring,"_", tk+offset, " ~ ",
                                               get_rhs(concurrentY=concurrentY, timepoint = tk#+1
                                                       , name_baseline_covariates = name_baseline_covariates,
                                                       name_time_covariates = name_time_covariates,
                                                       name_regimen = name_regimen, regimen = TRUE,
                                                       Markov = Markov, constant_variables = constant_variables))}
          else {})
      }
    })))
  }
  
  
  #NOTE: Below is different for simulation
  ## Note that A_k ~ V + L_0 + ... + L_(k-1) + A_(k-1) for k = 1,..., time_horizon, but A_0 ~ V + L_0
  ## i.e., regimen at baseline depends on additional baseline covariates, whereas in general, regimen depends
  ## on the previously observed covariates and regimen.
  ## The reason for this is we do not want to mistakenly assume that L_1 -> A_1 when in reality A_1 happens before L_1
  Qform <- unlist(lapply(1:time_horizon,function(tk){
    c(if(tk==1 & firstL){    paste0("Q.kplus1 ~ ", get_rhs(concurrentY=FALSE, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                                           name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                                           Markov = Markov, constant_variables = constant_variables))}
      else{},
      paste0("Q.kplus1 ~ ", get_rhs(concurrentY=concurrentY, timepoint = tk, name_baseline_covariates = name_baseline_covariates,
                                    name_time_covariates  = name_time_covariates, name_regimen = name_regimen, regimen = TRUE,
                                    Markov = Markov, constant_variables = constant_variables)))
  }))
  
  if(concurrentY){
    for(i in seq_along(gform)){
      if(any(RhsVars(gform[[i]]) %in% LhsVars(gform[[i]]))){
        gform[[i]]=gsub(paste0(" \\+ ",RhsVars(gform[[i]])[RhsVars(gform[[i]]) %in% LhsVars(gform[[i]])]),"", gform[[i]])
      }
    }
  }
  
  if(firstL){
    names(Qform)=c(paste0(name_time_covariates[1],"_1"),paste0(name_outcome,"_",1:time_horizon))
  }else{
    names(Qform)=paste0(name_outcome,"_",1:time_horizon)
  }
  list(gform = gform, Qform = Qform)
}  
## Adjust data to fit ltmle constraints
get_ltmle_data <- function(work_data, time_horizon,
                           name_outcome,
                           name_baseline_covariates,
                           name_time_covariates,
                           name_regimen,
                           name_censoring,
                           censored_label,
                           name_comp.event){
  time_grid = 0:time_horizon
  K = length(time_grid)
  
  if(length(name_censoring)>0){
    for(col in sapply(time_grid[-1], function(timepoint){paste0(name_censoring,"_",timepoint)})){
      set(work_data, j = col, value=ifelse(work_data[[col]]==censored_label,"censored","uncensored"))
      set(work_data, j = col, value=as.factor(work_data[[col]]))}
  }
  
  ## Manipulation of the event nodes
  if(length(name_regimen)==2){
    A_nodes = c(t(matrix(c(paste0(name_regimen[[1]],"_",time_grid[-K]),paste0(name_regimen[[2]],"_",time_grid[-K])), nrow = K-1)))
  } else {
    A_nodes = paste0(name_regimen[[1]],"_",time_grid[-K])
  }
  Y_nodes = paste0(name_outcome,"_",time_grid[-1])
  D_nodes = paste0(name_comp.event,"_",time_grid[-c(1,K)])
  C_nodes = paste0(name_censoring,"_",time_grid[-1])
  A_nodes_position = match(A_nodes,names(work_data))
  Y_nodes_position = match(Y_nodes,names(work_data))
  D_nodes_position = match(D_nodes,names(work_data))
  C_nodes_position = match(C_nodes,names(work_data))
  
  
  ## Adjust data depending on censoring/event/competing event with NA
  
  for(q in 1:(K-1)){
    if(q<(K-1)){
      has_outcome_or_death_and_censored = (((work_data[[Y_nodes_position[[q]]]]%in%"1")|(work_data[[D_nodes_position[[q]]]]%in%"1"))&
                                             (work_data[[C_nodes_position[[q]]]]%in%"censored"))
    } else{
      has_outcome_or_death_and_censored = ((work_data[[Y_nodes_position[[q]]]]%in%1)&(work_data[[C_nodes_position[[q]]]]%in%"censored"))
    }
    if(any(has_outcome_or_death_and_censored)){
      set(work_data,j=C_nodes_position[[q]],i=which(has_outcome_or_death_and_censored),value="uncensored")
    }
  }
  
  ## All nodes (except outcome and competing risk) should be NA after an event (outcome or death)
  
  if(time_horizon!= 1){
    for(k in Y_nodes_position[-(K-1)]){
      later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
      later_Y_nodes=intersect((k+1):Y_nodes_position[length(Y_nodes_position)],Y_nodes_position)
      if(any(has_outcome <- (work_data[[k]]%in%1))){
        for(l in later_nodes) {set(work_data,j=l,i=which(has_outcome),value=NA)}
        for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_outcome),value=1)}
      }
    }
    if(length(name_comp.event)>0){
      for(k in D_nodes_position){
        later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
        # Later outcome event nodes are set to 0
        later_Y_nodes=intersect((k+1):NCOL(work_data),Y_nodes_position)
        if(any(has_died <- (work_data[[k]]%in%1))){
          for(l in later_nodes) {set(work_data,j=l,i=which(has_died),value=NA)}
          for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_died),value=0)}
        }
      }
    }
    ## All nodes should be NA as soon as censoring has occurred
    if(length(name_censoring)>0){
      for(k in C_nodes_position){
        later_nodes=(k+1):NCOL(work_data)
        if(any(has_censored <- (work_data[[k]]%in%"censored"))){
          for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
        }
      }
    }
  }else{
    # 
    # time_horizon = 1 we set the outcome to NA in case of censored
    #
    if(length(name_censoring)>0){
      for(k in C_nodes_position){
        later_nodes=(k+1):NCOL(work_data)
        if(any(has_censored <- (work_data[[k]]%in%"censored"))){
          for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
        }
      }
    }
  }
  # Data at risk
  at.risk = list()
  at.risk[[1]] = work_data[,.N]
  if(time_horizon > 1){
    for(i in 2:time_horizon){
      if(length(name_censoring)>0 & length(name_comp.event)>0){
        at.risk[[i]] = sum((work_data[[Y_nodes_position[[i-1]]]]%in%0)
                           &(work_data[[C_nodes_position[[i-1]]]]%in%"uncensored")
                           &(work_data[[D_nodes_position[[i-1]]]]%in%0))
      } else {
        if(length(name_censoring)>0){
          at.risk[[i]] = sum((work_data[[Y_nodes_position[[i-1]]]]%in%0)
                             &(work_data[[C_nodes_position[[i-1]]]]%in%"uncensored"))
        } else{
          if(length(name_comp.event)>0){
            at.risk[[i]] = sum((work_data[[Y_nodes_position[[i-1]]]]%in%0)
                               &(work_data[[D_nodes_position[[i-1]]]]%in%0))
          } else{
            at.risk[[i]] = sum(work_data[[Y_nodes_position[[i-1]]]]%in%0)
          }
        }
      }
    }
    names(at.risk) = paste0("Number of persons at risk at time point ", c(2:time_horizon))
  }
  
  # L_nodes = c(name_baseline_covariates, sapply(time_grid, function(k) {paste0(c(name_time_covariates, name_comp.event), "_", k)}))
  L_nodes = c(sapply(time_grid, function(k) {paste0(c(name_comp.event,name_time_covariates), "_", k)}))
  L_nodes = L_nodes[match(L_nodes, names(work_data),nomatch = 0)!=0]
  
  if(length(name_censoring)==0) {C_nodes = NULL}
  
  list(data = work_data[],
       Anodes = A_nodes,
       Cnodes = C_nodes,
       Lnodes = L_nodes, 
       Ynodes = Y_nodes,
       at.risk = at.risk)
}
## rhs of formulas used for GLM, g-formulas and Q-formulas
get_rhs <- function(timepoint, name_baseline_covariates, name_time_covariates, name_regimen, regimen = TRUE, Markov = NULL,
                    constant_variables = NULL){
  form = paste(setdiff(name_baseline_covariates,constant_variables), collapse = "+")
  
  # name_time_covariates = paste0(NULL,if(timepoint!=0){name_time_covariates}else{}) # Include if A_0 ~ V and C_1 ~ V only
  
  if(length(name_time_covariates[name_time_covariates%in%Markov])>0) {
    form = paste(form, "+", paste(setdiff(sapply(name_time_covariates[name_time_covariates%in%Markov],
                                                 function(ntc) {paste0(ntc, "_", max(0, (timepoint - 1)))}),
                                          constant_variables),
                                  collapse = " + "))
  }
  if(length(name_time_covariates[!(name_time_covariates%in%Markov)])>0){
    form = paste(form, "+", paste(setdiff(sapply(name_time_covariates[!(name_time_covariates%in%Markov)],
                                                 function(ntc) {paste0(ntc, "_", 0:max(0, (timepoint - 1)))}),
                                          constant_variables),
                                  collapse = " + "))
  }
  if(regimen == TRUE) {
    form = paste(form, "+", paste(setdiff(sapply(name_regimen, function(nt) {paste0(nt, "_", 0:max(0, (timepoint - 1)))}),
                                          constant_variables),
                                  collapse = " + "))
  }
  form[]
}

## Subset rows and columns

get_subset_data <- function(work_data,
                            time_horizon,
                            subset_id = NULL,
                            subset_label = NULL,
                            name_baseline_covariates){
  
  time_grid = 0:time_horizon
  K = length(time_grid)
  if(length(subset_id)>0){
    work_data = work_data[pnr%in%subset_id]
  }
  # label the variables that are constant in the subset data
  same = sapply(work_data, function(x){length(unique(x))==1})
  if(sum(same)>0)
    constant_variables <- names(work_data)[same]
  else
    constant_variables <- NULL
  list(data = work_data[],
       subset_label = subset_label,
       name_baseline_covariates = name_baseline_covariates[name_baseline_covariates%in%names(work_data)],
       constant_variables= constant_variables)
}

ltmle.glmnet <- function(Y,
                         X,
                         newX,
                         family,
                         obsWeights,
                         id,
                         alpha = 1,
                         nfolds = 10,
                         selector="undersmooth",
                         nlambda = 100,
                         useMin = TRUE,
                         loss = "deviance",
                         ...){
  
  requireNamespace("glmnet")
  
  Xnames=attr(newX,"Xnames")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  FAM <- ifelse(length(unique(Y))>2,"gaussian","binomial")
  if (length(selector)>0&&selector=="undersmooth")
    uoh <- try(fit <- glmnet::glmnet(X,Y,weights = obsWeights,lambda = NULL,alpha = alpha,nlambda = nlambda,trace.it = 0L,family=FAM,...))
  else{
    # make sure that 
    if (any(duplicated(id))){
      id_data=data.table(id=id)
      foldid_data=data.table(id=unique(id),foldid=sample(1:nfolds,size=length(unique(id)),replace=TRUE))
      foldid=foldid_data[id_data,on="id"]
    }else{
      foldid=NULL
    }
    
    uoh <- try(fit <- glmnet::cv.glmnet(x = X,
                                        y = Y,
                                        weights = obsWeights,
                                        lambda = NULL,
                                        type.measure = loss,
                                        nfolds = nfolds,
                                        foldid=foldid,
                                        family = FAM,
                                        alpha = alpha,
                                        nlambda = nlambda,
                                        ...))
  }
  if (inherits(uoh,"try-error")) #browser()
    stop("ltmle.glmnet could not fit")
  if (length(selector)>0&&selector=="undersmooth"){
    selected.lambda <- fit$lambda[length(fit$lambda)]
    pred <- c(predict(fit, newx = newX, type = "response", s = fit$lambda[length(fit$lambda)]))
  } else{
    ## pred <- c(predict(fit$glmnet.fit, newx = newX, type = "response", s = ifelse(useMin,"lambda.min","lambda.1se")))
    if (useMin)
      selected.lambda <- fit$lambda.1se
    else
      selected.lambda <- fit$lambda.min
    pred <- c(predict(fit$glmnet.fit, newx = newX, type = "response", s = selected.lambda))
  }
  class(fit) <- c("ltmle.glmnet")
  if (length(selector)>0&&selector=="undersmooth"){
    beta <- data.table::data.table(beta=fit$beta[,NCOL(fit$beta),drop=TRUE])
  } else{
    bmat <- fit$glmnet.fit$beta
    beta <- data.table::data.table(beta=bmat[,match(fit$lambda.min,fit$lambda)])
  }
  if (length(Xnames)==NROW(beta)) beta=cbind(X=Xnames,beta)
  attr(fit,"selector") <- selector
  fit$selected_beta <- beta
  fit$selected.lambda <- selected.lambda
  out <- list(predicted.values = pred, fit = fit)
  return(out)
}

predict.ltmle.glmnet <- function(object,newX,...){
  selector <- attr(object,"selector")
  if (!is.matrix(newX)){
    newX <- model.matrix(~-1 + ., newX)
  }
  if (length(selector)>0&&selector=="undersmooth"){
    class(object)="glmnet"
    uha <- try(pred <- c(predict(object,
                                 newx = newX,
                                 type = "response",
                                 s = object$lambda[length(object$lambda)],
                                 ...)))
    if (inherits(uha,"try-error"))
      stop("Prediction of glmnet object failed in ltmle.glmnet.")
    else
      pred
  } else{
    pred <- c(predict(object$glmnet.fit,
                      newx = newX,
                      type = "response",
                      s = object$selected.lambda,
                      ... ))
  }
}
Ltmle <- function (data, Anodes, Cnodes = NULL, Lnodes = NULL, Ynodes,
                   survivalOutcome = NULL, Qform = NULL, gform = NULL, abar,
                   rule = NULL, gbounds = c(0.01, 1), Yrange = NULL, deterministic.g.function = NULL,
                   stratify = FALSE, SL.library = "glm", SL.cvControl = list(),
                   estimate.time = TRUE, gcomp = FALSE, iptw.only = FALSE, deterministic.Q.function = NULL,
                   variance.method = "tmle", observation.weights = NULL, id = NULL,info = NULL,verbose=FALSE)
{
  nn=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
  require(matrixStats)
  data <- CheckData(data)
  msm.inputs <- GetMSMInputsForLtmle(data, abar, rule, gform)
  inputs <- CreateInputs(data = data, Anodes = Anodes, Cnodes = Cnodes,
                         Lnodes = Lnodes, Ynodes = Ynodes, survivalOutcome = survivalOutcome,
                         Qform = Qform, gform = msm.inputs$gform, Yrange = Yrange,
                         gbounds = gbounds, deterministic.g.function = deterministic.g.function,
                         SL.library = SL.library, SL.cvControl = SL.cvControl,
                         regimes = msm.inputs$regimes, working.msm = msm.inputs$working.msm,
                         summary.measures = msm.inputs$summary.measures, final.Ynodes = msm.inputs$final.Ynodes,
                         stratify = stratify, msm.weights = msm.inputs$msm.weights,
                         estimate.time = estimate.time, gcomp = gcomp, iptw.only = iptw.only,
                         deterministic.Q.function = deterministic.Q.function,
                         variance.method = variance.method, observation.weights = observation.weights,
                         id = id, verbose = verbose)
  result <- LtmleFromInputs(inputs)
  result$call <- match.call()
  result$info <- result$call$info
  class(result) <- "Ltmle"
  return(result)
}
merge_and_sort_data <- function(time_horizon,
                                regimen_data,
                                outcome_data, 
                                baseline_data,
                                timevar_data, 
                                name_outcome,                    
                                name_regimen,
                                name_censoring = NULL,
                                censored_label = "censored",
                                name_comp.event = NULL,
                                test=FALSE){
  time_grid = 0:time_horizon
  K = length(time_grid)
  # the regimen may have two components (A and B) both are
  # identified based on the first (A)
  # here we select the specified element of a list of regimens
  setkey(regimen_data,pnr)
  outcome_data <- outcome_data[[name_outcome]]
  setkey(outcome_data,pnr)
  wide_data=outcome_data[regimen_data]
  #
  # deal with outcome/death/censored at index
  #
  # D_0 = match(paste0(name_comp.event,"_",0),names(wide_data))
  # C_0 = match(paste0(name_censoring,"_",0),names(wide_data))
  # wide_data = wide_data[!(wide_data[[D_0]]%in%1)&!(wide_data[[C_0]]%in%censored_label)]
  #
  # adding the baseline covariates
  #
  wide_data=baseline_data[wide_data, on = c("pnr")]
  # subset and sort data
  if (test>0){
    if (test==1||test==TRUE) test=min(20000,NROW(wide_data))
    message("\nSubsetting to ",test," randomly selected rows of the data.")
    work_data <- wide_data[sample(1:.N,size=test,replace=FALSE)]
  } else{
    work_data <- wide_data
  }
  # add time covariates
  # first remove outcome if overlap
  if (length((outcome_overlap <- grep(paste0(name_outcome,"_"),names(timevar_data)))
             >0))
    timevar_data <- timevar_data[,-outcome_overlap, with=FALSE]
  setkey(timevar_data,pnr)
  work_data=timevar_data[work_data, on = c("pnr")]
  
  name_time_covariates = unlist(lapply(grep("_0",names(timevar_data),value=TRUE),function(x){substring(x,0,nchar(x)-2)}))
  name_baseline_covariates = setdiff(names(baseline_data),"pnr")
  
  # sorting the variables for LTMLE
  work_data = work_data[,c("pnr", name_baseline_covariates,unlist(sapply(time_grid, function(timepoint){
    if(timepoint == 0){
      paste0(c(name_time_covariates, name_regimen),"_",timepoint)
    } else{
      if(timepoint != time_grid[K]){
        paste0(c(name_censoring, name_outcome, name_comp.event, name_time_covariates, name_regimen),"_",timepoint)
      } else {
        paste0(c(name_censoring, name_outcome),"_",timepoint)
      }
    }
  }))), with = FALSE]
  
  list(data = work_data[],
       name_baseline_covariates = name_baseline_covariates,
       name_time_covariates = name_time_covariates,
       name_regimen = name_regimen)
}

prepare_Ltmle <- function(regimen_data,
                          outcome_data,
                          baseline_data,
                          timevar_data,
                          time_horizon,
                          subset_id = NULL,
                          subset_label = NULL,
                          name_outcome,
                          name_regimen,
                          name_censoring = "Censored",
                          censored_label = "censored",
                          name_comp.event = "Dead",
                          Markov = NULL,
                          abar,
                          deterministic.Q.function = NULL,
                          SL.library,
                          test = FALSE) {
  ## Merge all data and order in correct order
  merged_data = merge_and_sort_data(time_horizon = time_horizon,
                                    regimen_data = regimen_data,
                                    outcome_data = outcome_data,
                                    baseline_data = baseline_data,
                                    timevar_data = timevar_data,
                                    name_outcome = name_outcome,
                                    name_regimen = name_regimen,
                                    name_censoring = name_censoring,
                                    censored_label = censored_label,
                                    name_comp.event = name_comp.event,
                                    test = test)
  
  ## Subsetting the data; This returns data in correct order according to time and without constant nodes
  subset_data = get_subset_data(work_data = merged_data$data,
                                time_horizon = time_horizon,
                                subset_id = subset_id,
                                subset_label = subset_label,
                                name_baseline_covariates = merged_data$name_baseline_covariates)
  
  ## Change data to fit into ltmle constraints; Censored should be factor with levels "uncensored" and "censored",
  ## all nodes occurring after censoring should be NA, all nodes (except outcome) occurring after an event (outcome or competing) should be NA
  ltmle_data = get_ltmle_data(subset_data$data,
                              time_horizon = time_horizon,
                              name_outcome = name_outcome,
                              name_baseline_covariates = subset_data$name_baseline_covariates,
                              name_time_covariates=merged_data$name_time_covariates,
                              name_regimen=name_regimen,
                              name_censoring=name_censoring,
                              censored_label=censored_label,
                              name_comp.event=name_comp.event)
  
  formulas = get_formulas(time_horizon = time_horizon,
                          work_data = ltmle_data$data,
                          name_outcome = name_outcome,
                          name_baseline_covariates = subset_data$name_baseline_covariates,
                          name_time_covariates = merged_data$name_time_covariates,
                          name_regimen = name_regimen,
                          name_censoring = name_censoring,
                          name_comp.event = name_comp.event,
                          Markov = Markov,
                          constant_variables = subset_data$constant_variables)
  
  ## abar
  if (missing(abar)){
    if(length(name_regimen)==2) {
      abar = list(rep(1:0,time_horizon),rep(0:1,time_horizon))}
    else
      abar <- list(rep(1,time_horizon), rep(0,time_horizon))
  }
  ## Deterministic Q function -- creates function indicating that competing risk event means that no event can happen
  dq <- function(data, current.node, nodes, called.from.estimate.g){
    death.index <- grep(paste0(name_comp.event, "_"),names(data))
    if(length(death.index)==0)stop("No death/terminal event node found")
    hist.death.index <- death.index[death.index < current.node]
    if(length(hist.death.index)==0)
      return(NULL)
    else{
      is.deterministic <- Reduce("+",lapply(data[,hist.death.index,drop=FALSE],
                                            function(dd){x=dd;x[is.na(dd)] <- 0;x}))>=1
      # should be unnecessary to exclude those who readily
      # have a missing value for death, but it does not hurt either
      is.deterministic[is.na(is.deterministic)] <- FALSE
      list(is.deterministic=is.deterministic, Q.value=0)
    }
  }
  
  ## Make sure that competing risk and deterministic.Q.function agree
  det.Q.function = NULL
  if(length(deterministic.Q.function)>0){
    det.Q.function = deterministic.Q.function
  }
  if(length(deterministic.Q.function)==0&length(name_comp.event)>0){
    if(length(grep(name_comp.event, names(ltmle_data$data), value = TRUE))>0){
      det.Q.function = dq
    }
  }
  
  ## Message about the time interval
  time_interval = NULL
  list(data = ltmle_data$data[],
       Qform = formulas$Qform,
       gform = formulas$gform,
       estimate.time = FALSE,
       Anodes = ltmle_data$Anodes,
       Cnodes = ltmle_data$Cnodes,
       Lnodes = ltmle_data$Lnodes,
       Ynodes = ltmle_data$Ynodes,
       survivalOutcome = TRUE,
       abar = abar,
       deterministic.Q.function = det.Q.function,
       SL.library = SL.library,
       info = list(outcome = name_outcome,
                   regimen = name_regimen,
                   baseline = subset_data$name_baseline_covariates,
                   timevar = merged_data$name_time_covariates,
                   subset_label = subset_data$subset_label,
                   order = merged_data$order,
                   time_horizon = time_horizon,
                   time_interval = time_interval,
                   at.risk = ltmle_data$at.risk,
                   constant_variables= subset_data$constant_variables,
                   Markov = Markov)
       ## formulas = paste("The Anode at time 0 depends on Lnodes at time 0,",
       ## "whereas the Anode at time k depends on Lnodes and Anodes up to time k-1"))
  )
}
print.Ltmle <- function(x,...){
  cat("\nFitted object obtained with augmented Ltmle.
Use names() and summary() to see results")
}
### print.prepare_Ltmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 19 2023 (15:57) 
## Version: 
## Last-Updated: May 22 2023 (09:25) 
##           By: Thomas Alexander Gerds
##     Update #: 61
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
print.prepare_Ltmle <- function(x,...){
  cat("\nPrepared data for analysis with Ltmle.\n")
  order = ifelse(x$info$order_YC,"The events variables are ordered by: outcome, censoring",
                 "The event variables are ordered by: censoring, outcome")
  cat(order,"\n")
  cat("Sum of treated and number of events by time: \n\n")
  Time = x$info$time_grid
  Treat = c(sapply(x$Anodes,function(a){sum(x$data[[a]],na.rm = TRUE)}),NA)
  Outcome = c(0,sapply(x$Ynodes,function(y){sum(x$data[[y]],na.rm = TRUE)}))
  if (length(x$Cnodes)>0){
    Cens = c(0,sapply(x$Cnodes,function(c){sum(x$data[[c]] == "censored",na.rm = TRUE)}))
    Cens[length(Time)] = NA
  }
  if (length(x$info$Dnodes)>0){
    CR = c(0,sapply(x$info$Dnodes,function(d){sum(x$data[[d]],na.rm = TRUE)}),NA)
  }
  N.G = N.Q = rep(NROW(x$data),length(Time))
  if (length(Time)>1){
    N.G <- N.G-c(0,sapply(2:(length(Time) -1),function(j){
      sum((x$data[[x$Ynodes[j-1]]] == 1
           | (x$data[[x$Cnodes[j-1]]] == "censored")|is.na(x$data[[x$Cnodes[j-1]]])
           | (x$data[[x$info$Dnodes[j-1]]] == 1)),na.rm = TRUE)
    }),NA)
    if (x$info$order_YC){
      N.Q <- N.Q-c(0,sapply(2:(length(Time)-1),function(j){
        sum(((x$data[[x$Ynodes[j-1]]] == 1)
             | (x$data[[x$Cnodes[j-1]]] == "censored")|is.na(x$data[[x$Cnodes[j-1]]])
             | (x$data[[x$info$Dnodes[j-1]]] == 1)),na.rm = TRUE)
      }),NA)
    }else{
      N.Q <- N.Q-c(sum(x$data[[x$Cnodes[1]]] == "censored",na.rm = TRUE),sapply(2:(length(Time)-1),function(j){
        sum(((x$data[[x$Cnodes[j]]] == "censored")|is.na(x$data[[x$Cnodes[j]]])
             | (x$data[[x$Ynodes[j-1]]] == 1) 
             | (x$data[[x$info$Dnodes[j-1]]] == 1)),na.rm = TRUE)
      }),NA)
    }
    N.Q[length(Time)] = N.Q[length(Time)-1]-Outcome[length(Time)]
  }
  X = data.table(Time = Time,
                 NumberAtrisk.G = N.G,
                 NumberAtrisk.Q = N.Q,
                 Treat = Treat,
                 Outcome = Outcome)
  if (length(x$info$Dnodes)>0){X[,CompetingEvents := CR]}
  if (length(x$Cnodes)>0){X[,Censored := Cens]}
  print(X[])
  invisible(X[])
}

######################################################################
### print.prepare_Ltmle.R ends here
######################################################################
summary.Ltmle <- function(object,estimator,...){
    nn=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
    if (missing(estimator))
        if (object$gcomp) estimator = "gcomp" else estimator = "tmle"
    summary_ltmle <- function (object, estimator = ifelse(object$gcomp, "gcomp", "tmle"),
                               ...)
    {
        IC.variance <- var(object$IC[[estimator]])
        if (estimator == "tmle" && !is.null(object$variance.estimate)) {
            v <- max(IC.variance, object$variance.estimate)
        }
        else {
            v <- IC.variance
        }
        variance.estimate.ratio = v/IC.variance
        if (object$binaryOutcome) {
            CIBounds <- c(0, 1)
        }
        else {
            CIBounds <- c(-Inf, Inf)
        }
        treatment <- GetSummary(list(long.name = NULL, est = object$estimates[estimator],
                                     gradient = 1, log.std.err = FALSE, CIBounds = CIBounds),
                                v,
                                n = length(object$IC[[estimator]]))
        ans <- list(treatment = treatment, call = object$call, estimator = estimator,
                    variance.estimate.ratio = variance.estimate.ratio)
        class(ans) <- "summary.ltmle"
        return(ans)
    }
    summi <- function(x,target){
        with(x,data.table::data.table(
                               Target_parameter=target,
                               Estimator=estimator,
                               estimate=estimate,
                               std.err=std.dev,
                               lower=CI[[1]],
                               upper=CI[[2]],
                               pvalue=pvalue))
    }
    if (length(object$estimates)>0){
        x=summary_ltmle(object,estimator=estimator)
        risk = summi(x=x$treatment,"Risk(A=1)")
        risk
    }else{
        x= summary.ltmleEffectMeasures(object,estimator=estimator)
        treatment = summi(x=x$effect.measures$treatment,"Risk(A=1)")
        control = summi(x=x$effect.measures$control,"Risk(A=0)")
        ate = summi(x=x$effect.measures$ATE,"ATE")
        RR = summi(x=x$effect.measures$RR,"RelativeRisk")
        out=rbind(treatment,control,ate,RR)
        out
    }
}


BinaryToCensoring <-
  function (is.censored, is.uncensored) 
  {
    if (!xor(missing(is.censored), missing(is.uncensored))) 
      stop("exactly one of is.censored and is.uncensored must be passed")
    calling.name <- names(sys.call(0))[2]
    if (length(calling.name) == 0 || !calling.name %in% c("is.censored", 
                                                          "is.uncensored")) 
      stop("the argument to BinaryToCensoring must be completely named - see ?BinaryToCensoring")
    if (missing(is.uncensored)) {
      is.uncensored <- !is.censored
    }
    if (!all(is.uncensored %in% c(0, 1, NA))) 
      stop("the argument to BinaryToCensoring should be binary (0, 1, or NA) or logical")
    y <- character(length(is.uncensored))
    y[is.uncensored == 0] <- "censored"
    y[is.uncensored == 1] <- "uncensored"
    y[is.na(is.uncensored)] <- NA
    return(factor(y))
  }



d <- structure(list(sex = structure(c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 
                                      2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 
                                      2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 
                                      2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 
                                      1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 
                                      2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 
                                      1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 
                                      1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 
                                      1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 
                                      1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                      1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 
                                      2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                      1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 
                                      2L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 
                                      2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 
                                      2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 
                                      1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 
                                      1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 
                                      2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 
                                      2L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                      2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 
                                      1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 
                                      1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 
                                      1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
                                      2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 
                                      1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 
                                      1L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                      1L, 2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 
                                      1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 
                                      1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 
                                      2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 
                                      1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 
                                      2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 
                                      2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 
                                      2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
                                      2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                      2L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 
                                      2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 
                                      2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 
                                      1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 
                                      2L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                      2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 
                                      2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 
                                      2L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 
                                      1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 
                                      1L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 
                                      2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 
                                      2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 
                                      2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 2L, 
                                      2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 
                                      1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 
                                      2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 
                                      2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                      2L, 2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 
                                      2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 
                                      2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 
                                      2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 
                                      2L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 
                                      2L), levels = c("Female", "Male"), class = "factor"), hypertension_0 = c(0, 
                                                                                                               0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 
                                                                                                               0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                               1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 
                                                                                                               0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                               0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                               1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 
                                                                                                               1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 
                                                                                                               0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                               1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 
                                                                                                               0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 
                                                                                                               0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
                                                                                                               0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 
                                                                                                               0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 
                                                                                                               0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0), GLP1RA_0 = c(0, 0, 1, 1, 
                                                                                                                                                                 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 
                                                                                                                                                                 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 
                                                                                                                                                                 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 
                                                                                                                                                                 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 
                                                                                                                                                                 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0), dementia_1 = c(0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                            0, 0, 0, 0, 0, 0, 0), hypertension_1 = c(0, 0, 1, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 
                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 
                                                                                                                                                                                                                                                     1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                     0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 
                                                                                                                                                                                                                                                     0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                     1, 1, 1, 1, 0, 0), GLP1RA_1 = c(0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 
                                                                                                                                                                                                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 
                                                                                                                                                                                                                                                                                     0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                     0, 0, 0), dementia_2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                              0), Censored_1 = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                           2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                           1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
                                                                                                                                                                                                                                                                                                              ), levels = c("uncensored", "censored"), class = "factor"), Censored_2 = structure(c(1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L), levels = c("uncensored", "censored"
                                                                                                                                                                                                                                                                                                                                                                                                   ), class = "factor"), Dead_1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0), Dead_2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     -1000L), class = c("data.table", "data.frame"))