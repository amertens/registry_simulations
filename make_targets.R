


library(here)
library(targets)


source(paste0(here::here(),"/_targets.R"))
tar_make()




tar_visnetwork(targets_only = TRUE)


targets::tar_meta(fields = warnings, complete_only = TRUE)


#check
d <- tar_read(sim_data)
d

