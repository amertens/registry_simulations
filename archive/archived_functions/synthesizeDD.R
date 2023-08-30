### synthesizeDD.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  3 2021 (15:20) 
## Version: 
## Last-Updated: Nov  5 2021 (17:47) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Synthesizing longitudinal diabetes dementia followup data
##'
##' A sequence of logistic regression models Danmark Statistics
##' @title Synthesizing longitudinal diabetes dementia followup data
##' @param coefficients Intercepts and regression coefficients (log-odds-ratios) 
##' @return \code{lvm} object for simulation
##' @seealso \code{lvm}, \code{distribution}, \code{regression}, \code{sim}
##' @examples
##' library(lava)
##' library(data.table)
##' cc <- fread("data/coefficients.txt")
##' u <- synthesizeDD(cc)
##' d <- sim(u,1000)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>

synthesizeDD <- function(coefficients, A_name = "glp1", A=NULL){
  requireNamespace("lava")
  coefficients <- data.table(coefficients)
  
  if(!is.null(A)){
    if(A==1){ coefficients$`(Intercept)`[grepl(A_name, coefficients$var)] <- 99999999999}
    if(A==0){ coefficients$`(Intercept)`[grepl(A_name, coefficients$var)] <- -99999999999}
  }
  
  
  XNAMES <- names(coefficients)[-(1:3)]
  BETA <- coefficients[,-(1:3),with=0L]
  INTERCEPT <- coefficients[["(Intercept)"]]
  # empty lava model for simulation
  m <- lvm()
  distribution(m,"age_base") <- normal.lvm(mean=70,sd=10)
  distribution(m,"sex") <- binomial.lvm(p=0.4)
  m <- addvar(m,"ie_type")
  m <- addvar(m,"code5txt")
  m <- addvar(m,"quartile_income")
  # loop across time and variables
  for (j in 1:NROW(coefficients)){
    V <- coefficients$var[j]
    beta <- unlist(BETA[j,])
    X <- XNAMES[!is.na(beta)]
    beta <- beta[!is.na(beta)]
    # add V ~ Intercept + beta X
    distribution(m,V) <- binomial.lvm()
    intercept(m,V) <- INTERCEPT[j]
    regression(m,from=X,to=V) <- beta
  }
  class(m) <- c("synthesizeDD",class(m))
  m
}


