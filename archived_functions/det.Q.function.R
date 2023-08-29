det.Q.function<-function(data, current.node, nodes, called.from.estimate.g){
    compete.index <- grep("^Dead",names(data))
    hist.compete.index <- compete.index[compete.index < current.node]
    if(length(hist.compete.index)==0)
        return(NULL)
    else{
        is.deterministic <- Reduce("+",lapply(data[,hist.compete.index,drop=FALSE],
                                              function(dd){x=dd;x[is.na(dd)] <- 0;x}))>=1
        is.deterministic[is.na(is.deterministic)] <- FALSE
        list(is.deterministic=is.deterministic, Q.value=0)
    }
}
