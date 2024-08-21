

shh <-try(setwd("C:/Users/andre/Documents/jici/registry_simulations/"),silent = TRUE)
shh <-try(setwd("~/research/Methods/registry_simulations/"),silent = TRUE)

nn=lapply(list.files("./functions/", full.names = TRUE, recursive=TRUE), source)


#calculate truth
set.seed(12345)

truth=calc_realistic_truth(A_name = "glp1", nsamp=100000, return_data=FALSE)
truth[10,]

#need to run many times to get a good estimate of the truth and average here:
#average_truth()