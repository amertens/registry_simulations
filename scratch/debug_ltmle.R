
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

# function (data, Anodes, Cnodes = NULL, Lnodes = NULL, Ynodes, 
#           survivalOutcome = NULL, Qform = NULL, gform = NULL, abar, 
#           rule = NULL, gbounds = c(0.01, 1), Yrange = NULL, deterministic.g.function = NULL, 
#           stratify = FALSE, SL.library = "glm", SL.cvControl = list(), 
#           estimate.time = TRUE, gcomp = FALSE, iptw.only = FALSE, deterministic.Q.function = NULL, 
#           variance.method = "tmle", observation.weights = NULL, id = NULL) 
# {
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
                         estimate.time = estimate.time, gcomp = gcomp, iptw.only = F, 
                         deterministic.Q.function = deterministic.Q.function, 
                         variance.method = variance.method, observation.weights = observation.weights, 
                         id = id,verbose=FALSE)
  results <- ltmle:::LtmleFromInputs(inputs)
  results_Ltmle <- LtmleFromInputs(inputs)
  
  
  inputs2 <- CreateInputs(data = data, Anodes = Anodes, Cnodes = Cnodes, 
                         Lnodes = Lnodes, Ynodes = Ynodes, survivalOutcome = survivalOutcome, 
                         Qform = Qform, gform = msm.inputs$gform, Yrange = Yrange, 
                         gbounds = gbounds, deterministic.g.function = deterministic.g.function, 
                         SL.library = "SL.mean", SL.cvControl = SL.cvControl, 
                         regimes = msm.inputs$regimes, working.msm = msm.inputs$working.msm, 
                         summary.measures = msm.inputs$summary.measures, final.Ynodes = msm.inputs$final.Ynodes, 
                         stratify = stratify, msm.weights = msm.inputs$msm.weights, 
                         estimate.time = estimate.time, gcomp = gcomp, iptw.only = F, 
                         deterministic.Q.function = deterministic.Q.function, 
                         variance.method = variance.method, observation.weights = observation.weights, 
                         id = id,verbose=FALSE)
  results2 <- ltmle:::LtmleFromInputs(inputs2)
  results2_Ltmle <- LtmleFromInputs(inputs2)
  
  summary(results)$effect.measures$ATE$estimate
  summary(results2)$effect.measures$ATE$estimate
  summary(results_Ltmle)$effect.measures$ATE$estimate
  summary(results2_Ltmle)$effect.measures$ATE$estimate
  
  summary(results, estimator="iptw")$effect.measures$ATE$estimate
  summary(results2, estimator="iptw")$effect.measures$ATE$estimate
  summary(results_Ltmle, estimator="iptw")$effect.measures$ATE$estimate
  summary(results2_Ltmle, estimator="iptw")$effect.measures$ATE$estimate
  
  #the difference is the package versus loaded LtmleFromInputs- check into the differences
  
  msm.result <- ltmle:::MainCalcs(inputs)
  msm.result2 <- ltmle:::MainCalcs(inputs2)
  msm.result_Ltmle <- MainCalcs(inputs)
  msm.result2_Ltmle <- MainCalcs(inputs2)  

  msm.result$beta.iptw
  msm.result2$beta.iptw
  msm.result_Ltmle$beta.iptw
  msm.result2_Ltmle$beta.iptw
  
  msm.result$beta
  msm.result2$beta
  msm.result_Ltmle$beta
  msm.result2_Ltmle$beta  
  
  
  g.list <- EstimateG(inputs)
  g.list2 <- ltmle:::EstimateG(inputs2)
  g.list_Ltmle <- EstimateG(inputs)
  g.list2_Ltmle <- EstimateG(inputs2)
  
  g.list$fit[[2]]$GLP1RA_1[2,]
  g.list2$fit[[2]]$GLP1RA_1
  g.list_Ltmle$fit[[2]]$GLP1RA_1[2,]
  g.list2_Ltmle$fit[[2]]$GLP1RA_1

  all.msm.weights <- ltmle:::GetMsmWeights(inputs)
  all.msm.weights2 <- ltmle:::GetMsmWeights(inputs2)
  all.msm.weights_Ltmle <- GetMsmWeights(inputs)
  all.msm.weights2_Ltmle <- GetMsmWeights(inputs2)
  
  # iptw <- ltmle:::CalcIPTW(inputs, g.list$cum.g, all.msm.weights)
  # iptw2 <- ltmle:::CalcIPTW(inputs2, g.list$cum.g, all.msm.weights)
  # iptw_Ltmle <- CalcIPTW(inputs, g.list$cum.g, all.msm.weights)
  # iptw2_Ltmle <- CalcIPTW(inputs2, g.list$cum.g, all.msm.weights)
  # 
  # iptw$beta
  # iptw2$beta
  # iptw_Ltmle$beta
  # iptw2_Ltmle$beta
  
  #NOTE! The difference is how the CalcIPTW handles g weights
  
  iptw <- ltmle:::CalcIPTW(inputs, g.list$cum.g, all.msm.weights)
  iptw2 <- ltmle:::CalcIPTW(inputs2, g.list2$cum.g, all.msm.weights)
  iptw_Ltmle <- CalcIPTW(inputs, g.list_Ltmle$cum.g, all.msm.weights)
  iptw2_Ltmle <- CalcIPTW(inputs2, g.list2_Ltmle$cum.g, all.msm.weights)

  iptw$beta
  iptw2$beta
  iptw_Ltmle$beta
  iptw2_Ltmle$beta
  
  msm.weights=all.msm.weights
  cum.g=g.list$cum.g
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
  m.glm_Ltmle <- ltmle.glm(formula(inputs$working.msm), family = quasibinomial(), data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, center = FALSE)))
  m.glm2_Ltmle <- ltmle.glm(formula(inputs2$working.msm), family = quasibinomial(), data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, center = FALSE)))
  m.glm <- ltmle:::ltmle.glm(formula(inputs$working.msm), family = quasibinomial(), data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, center = FALSE)))
  m.glm2 <- ltmle:::ltmle.glm(formula(inputs2$working.msm), family = quasibinomial(), data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, center = FALSE)))
  
  summary(m.glm_Ltmle)
  summary(m.glm2_Ltmle)
  summary(m.glm)
  summary(m.glm2)
  
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  #ltmle:::CalcIPTW
  

    #note! this differs between ltmle and Ltmle
    #m.glm <- ltmle.glm(formula(inputs$working.msm), family = quasibinomial(), data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, center = FALSE)))
    
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
  
##XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  #CalcIPTW
    


    m.glm <- ltmle.glm(formula(inputs$working.msm), family = quasibinomial(), 
                       data = data.frame(Y = Y.vec, X.mat, weight.vec), weights = as.vector(scale(weight.vec, 
                                                                                                  center = FALSE)))

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

