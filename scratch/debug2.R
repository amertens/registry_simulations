
rm(list=ls())
gc()
shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)
nn=lapply(list.files("./Ltmle/Augmentation/", full.names = TRUE, recursive=TRUE), source)

load(here::here("my_environment.RData"))

pl

data = pl$data
Anodes = pl$Anodes
Cnodes = pl$Cnodes
Lnodes = pl$Lnodes
Ynodes = pl$Ynodes
survivalOutcome = pl$survivalOutcome
Qform = pl$Qform
gform = pl$gform
abar = pl$abar
rule = pl$rule
gbounds = pl$gbounds
Yrange = pl$Yrange
deterministic.g.function = NULL
stratify = FALSE
SL.library = "glm"
SL.cvControl = list()
estimate.time = TRUE
gcomp = FALSE
iptw.only = FALSE
deterministic.Q.function = pl$deterministic.Q.function
variance.method = "ic"
observation.weights = NULL
id = NULL
info = NULL
verbose=FALSE

  nn=lapply(list.files("./Ltmle/R/", full.names = TRUE, recursive=TRUE), source)
  require(matrixStats)
  data <- CheckData(data)
  msm.inputs <- GetMSMInputsForLtmle(data, abar, rule, gform)
  
  
  inputs <- CreateInputs(data = data, Anodes = Anodes, Cnodes = Cnodes,
                         Lnodes = Lnodes, Ynodes = Ynodes, survivalOutcome = survivalOutcome,
                         Qform = Qform, gform = msm.inputs$gform, Yrange = Yrange,
                         gbounds = gbounds, deterministic.g.function = deterministic.g.function,
                         SL.library = "glm", SL.cvControl = SL.cvControl,
                         regimes = msm.inputs$regimes, working.msm = msm.inputs$working.msm,
                         summary.measures = msm.inputs$summary.measures, final.Ynodes = msm.inputs$final.Ynodes,
                         stratify = stratify, msm.weights = msm.inputs$msm.weights,
                         estimate.time = FALSE, gcomp = gcomp, iptw.only = iptw.only,
                         deterministic.Q.function = deterministic.Q.function,
                         variance.method = variance.method, observation.weights = observation.weights,
                         id = id, verbose = verbose)
  
  inputs2 <- CreateInputs(data = data, Anodes = Anodes, Cnodes = Cnodes,
                         Lnodes = Lnodes, Ynodes = Ynodes, survivalOutcome = survivalOutcome,
                         Qform = Qform, gform = msm.inputs$gform, Yrange = Yrange,
                         gbounds = gbounds, deterministic.g.function = deterministic.g.function,
                         SL.library = "glmnet", SL.cvControl =  list(selector="undersmoothed",alpha=1),
                         regimes = msm.inputs$regimes, working.msm = msm.inputs$working.msm,
                         summary.measures = msm.inputs$summary.measures, final.Ynodes = msm.inputs$final.Ynodes,
                         stratify = stratify, msm.weights = msm.inputs$msm.weights,
                         estimate.time = FALSE, gcomp = gcomp, iptw.only = iptw.only,
                         deterministic.Q.function = deterministic.Q.function,
                         variance.method = variance.method, observation.weights = observation.weights,
                         id = id, verbose = verbose)
  #result <- LtmleFromInputs(inputs)
  
  msm.result <- LtmleMSMFromInputs(inputs)
  msm.result$fit$g[[1]]$GLP1RA_1

  msm.result2 <- LtmleMSMFromInputs(inputs2)
  msm.result2$fit$g[[1]]$GLP1RA_1
  
  msm.result$beta.iptw
  msm.result2$beta.iptw
  
  summary(msm.result$IC.iptw[, 1])
  summary(msm.result2$IC.iptw[, 1])
  
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
  }else{
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
  #}

  class(result) <- "Ltmle"
