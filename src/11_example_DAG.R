library(dagitty)
library(ggdag)
library(ggraph)
library(cowplot)
library(here)
library(tidyverse)

# Create DAG structure with simplified paths
dag <- dagitty('
  dag {
    L0 -> A1 -> A2
    L0 -> Y1 -> Y2 -> Y3
    L0 -> Y2
    L0 -> Y3
    L0 -> A1
    L0 -> A2
    A1 -> Y2
    A1 -> Y3
    A2 -> Y3
    Y1 -> A2
    Y1 -> Y3
    Y2 -> Y3
  }
')

# Set coordinates with improved spacing
coordinates(dag) <- list(
  x = c(L0 = 0,    # Initial node
        A1 = 1,    # First layer
        Y1 = 1, 
        A2 = 3,    # Second layer
        Y2 = 3,
        Y3 = 5),   # Final layer
  y = c(L0 = 2,    # Center initial node
        A1 = 2.6,  # Top track
        Y1 = 1.4,  # Bottom track
        A2 = 2.6,  # Top track
        Y2 = 1.4,  # Bottom track
        Y3 = 1.8)    # Center final node
)

# Create plot with improved aesthetics
p <- dag %>% 
  tidy_dagitty(layout = "auto", seed = 12345) %>%
  arrange(name) %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "gray50", 
                 edge_width = 0.5,
                 arrow_directed = grid::arrow(length = grid::unit(5, "pt"))) +
  geom_dag_text(
    colour = "black",
    parse = TRUE,
    label = c(
      expression(A[1]),
      expression(A[t-1]),
      expression(W),
      expression(L[1]),
      expression(L[t-1]),
      expression(L[t])
    ),
    size = 4
  ) +
  theme_dag_blank() +
  theme(plot.margin = margin(20, 20, 20, 20))

# Display plot
p

# Save the plot
ggsave(p, file=paste0(here(),"/figures/longitudinal_dag.png"), width = 5, height = 3)