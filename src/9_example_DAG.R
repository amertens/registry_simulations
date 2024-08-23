
library(dagitty)
library(ggdag)
library(ggraph)
library(cowplot)
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
    D1 -> D2
    A2 -> D2
    L1 -> D2
    D2 -> Y2
  }
')

# Set coordinates for better visualization
coordinates(dag) <- list(
  x=c(L0=0, A1=1, L1=1.5, A2=3.5, Y1=3, Y2=5, D1=3, D2=5),
  y=c(L0=0, A1=1, L1=0, A2=1, Y1=2, Y2=2, D1=-1, D2=-1)
)


# # Plot the DAG
# p <- ggdag(dag, node=T) +
#   theme_dag_blank() +
#   ggtitle("Simplified DAG for Longitudinal Analysis of GLP-1RA Effect on Dementia")
# p
# 
# temp<-dag %>% 
#   tidy_dagitty(layout = "auto", seed = 12345) %>%
#   arrange(name)
# 
# unique(temp$data$name)

p <- dag %>% 
  tidy_dagitty(layout = "auto", seed = 12345) %>%
  arrange(name) %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  #geom_dag_point() +
  geom_dag_edges() +
  geom_dag_text(colour="black", parse = TRUE, label = c(expression(tilde(A)(t)),
                                        expression(tilde(A)(t+1)),
                                        expression(D(t)),
                                        expression(D(t+1)),
                                        expression(tilde(L)(0)), 
                                        expression(tilde(L)(t)), 
                                        expression(Y(t)),
                                        expression(Y(t+1)))) +  
  theme_dag_blank()
p  

# Save the plot
ggsave(p, file=paste0(here(),"/figures/longitudinal_dag.png"), width = 5, height = 3)
