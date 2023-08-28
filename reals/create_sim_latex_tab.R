

create_sim_latex_tab <- function(res_table){
  
    # identify index of rows to highlight
    row.i.1 <- which(res_table$filenames=="sim_res_DetQ_ic_v3")
    
    res_table <- res_table %>% select(filenames,  iptw,
                                      estimator, Qint,
                                      bias,variance,bias_se_ratio,oracle.coverage) %>%
      rename(Algorithm=estimator,  `Q-int`=Qint,  Estimator=iptw, `Bias/SE`=bias_se_ratio,
             Bias=bias, Variance=variance, `Oracle coverage`=oracle.coverage)
    
    row.i.1 <- which(res_table$filenames=="sim_res_DetQ_ic_v3")
    
    
    #save for html file
    res_table_protective <- res_table
    res_table_protective_raw <- res_diff_raw

  res_table <- res_table %>% subset(., select = -c(filenames))
  
  print(as.data.frame(res_table))
  
  res_xtable <- res_table %>%
    knitr::kable(
      format = "latex",
      align = "l",
      booktabs = TRUE,
      longtable = TRUE,
      linesep = "",
      digits =6) %>%
    kableExtra::kable_styling(
      position = "left"#,
      # latex_options = c("striped", "repeat_header"),
      # stripe_color = "gray!15"
    ) %>%
    kable_styling()%>%
    row_spec(row.i.1-1, hline_after = T) %>%
    row_spec(row.i.1, bold=T,hline_after = T)
  
  
  #save_kable(res_xtable, file="C:/Users/andre/Documents/jici/diab-dementia-server-code/tables/sim_results_table.tex",float = FALSE)

}
