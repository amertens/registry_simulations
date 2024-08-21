
library(dagitty)
library(ggdag)
library(here)
library(tidyverse)


# Create the DAG

# Create the DAG
dag <- dagitty('
  dag {
    L0 -> A1 -> Y1 -> Y2
    L0 -> L1 -> Y1 -> Y2
    L0 -> L1 -> A2 -> Y2
    A1 -> L1
    L1 -> A2
    L0 -> A2
    A1 -> A2
    L0 -> Y1
    L0 -> Y2
    A1 -> Y2
    A1 -> D1
    L1 -> D1
    D1 -> Y1
    D1 -> Y2
    A2 -> D2
    L1 -> D2
    D2 -> Y2
  }
')

# Set coordinates for better visualization
coordinates(dag) <- list(
  x=c(L0=0, A1=1, L1=2, A2=3, Y1=2, Y2=4, D1=1, D2=3),
  y=c(L0=0, A1=1, L1=0, A2=1, Y1=2, Y2=2, D1=-1, D2=-1)
)


# Plot the DAG
p <- ggdag(dag) +
  theme_dag() +
  ggtitle("Simplified DAG for Longitudinal Analysis of GLP-1RA Effect on Dementia")
p

# Save the plot
ggsave(p, file=paste0(here(),"/figures/longitudinal_dag.png"), width = 10, height = 6)
