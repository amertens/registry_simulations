### simulate_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 17 2023 (13:50) 
## Version: 
## Last-Updated: Jul 17 2023 (14:14) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log: 
#----------------------------------------------------------------------
## 
### Code:
simulate_data <- function(lava_model,n){
    d = lava::sim(lava_model,n)
    data.table::setDT(d)
    setnames(d,"sexMale","sex")
    d[,sex:=factor(sex,levels=c("0","1"),labels=c("Female","Male"))]
    cnames = grep("Censored_",names(d),value = TRUE)
    for (c in cnames)
        set(d,j = c,value = factor(d[[c]],levels=c(1,0),labels=c("uncensored","censored")))
    # from integer to numeric
    isINT=sapply(names(d),function(n)is.integer(d[[n]]))
    d[,(names(d)[isINT]):=lapply(.SD, as.numeric), .SDcols = names(d)[isINT]]
    d[]
}

######################################################################
### simulate_data.R ends here
