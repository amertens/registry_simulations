
det.Q.function <- function(data, current.node, nodes, called.from.estimate.g){
  death.index <- grep("death_",names(data))
  if(length(death.index)==0)stop("no node found")
  hist.death.index <- death.index[death.index < current.node]
  if(length(hist.death.index)==0)return(NULL)
  if(length(hist.death.index)==1){
    is.deterministic <- data[,hist.death.index]==1
  } else {
    is.deterministic <- apply(data[,hist.death.index]==1,1,any)
  }#if death before, remove from fitting
  is.deterministic[is.na(is.deterministic)] <- F
  return(list(is.deterministic=is.deterministic, Q.value=0))
  
}