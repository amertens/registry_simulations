

clean_ltmle_res <- function(fit, analysis_name="", iteration=NULL){

  res = summary(fit)
  
  #add in
  res.iptw = summary(fit,"iptw")
  
  res=bind_rows(res, res.iptw)
  # res$Qform = list(fit$formulas$Qform)
  # res$gform = list(fit$formulas$gform)
  res$analysis_name=analysis_name
  res$iteration=iteration
  res=data.frame(res)
  return(res)  
}