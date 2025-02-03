library(ltmle)

# Helper functions
rexpit <- function(x) rbinom(n=length(x), size=1, prob=(x))

# Generate data for n subjects across 3 time points
n <- 1000
set.seed(123)

# Baseline covariates
W <- data.frame(
  W1 = rnorm(n),  # continuous covariate
  W2 = rbinom(n, 1, 0.4)  # binary covariate
)

# Generate time-varying data
data <- W

# For each time point
for(t in 1:3) {
  # Treatment
  if(t == 1) {
    prob_a <- plogis(-1 + 0.3*W$W1 + 0.5*W$W2)
  } else {
    prob_a <- plogis(-1 + 0.3*data[[paste0("A", t-1)]])
  }
  data[[paste0("A", t)]] <- rexpit(prob_a)
  
  # Death (competing risk)
  if(t == 1) {
    prob_d <- plogis(-6 + 0.2*W$W1 - 0.5*data[[paste0("A", t)]])
  } else {
    prob_d <- plogis(-6 + 0.2*data[[paste0("A", t)]] - 0.5*data[[paste0("A", t-1)]])
  }
  data[[paste0("D", t)]] <- rexpit(prob_d)
  
  # Outcome
  if(t == 1) {
    prob_y <- plogis(-5 + 0.3*W$W1 + data[[paste0("A", t)]])
  } else {
    prob_y <- plogis(-5 + 0.3*data[[paste0("A", t)]] + 0.2*data[[paste0("A", t-1)]])
  }
  data[[paste0("Y", t)]] <- rexpit(prob_y)
}

table(data$Y1, data$A1)
table(data$D1, data$A1)

table(data$Y1)
table(data$Y2)
table(data$Y3)
table(data$Y3, data$A3)

table(data$D1)
table(data$D2)
table(data$D3)
table(data$D3, data$A3)

table(data$A1)
table(data$A2)
table(data$A3)

# Process data for survival outcome format
# Function to properly set subsequent values to 1 after any event
process_survival_data <- function(data) {
  n <- nrow(data)
  Ynodes <- c("D1", "Y1", "D2", "Y2", "D3", "Y3")
  
  for(i in 1:n) {
    had_event <- FALSE
    for(node in Ynodes) {
      if(had_event) {
        data[i, node] <- 1
      } else if(data[i, node] == 1) {
        had_event <- TRUE
      }
    }
  }
  return(data)
}

# Apply the processing
data <- process_survival_data(data)



# Function to check data patterns
check_data_patterns <- function(data) {
  # Get indices of first events
  n <- nrow(data)
  Ynodes <- c("D1", "Y1", "D2", "Y2", "D3", "Y3")
  Anodes <- paste0("A", 1:3)
  
  # Pick a few rows where events occurred
  for(i in 1:min(5, n)) {
    # Find first event
    first_event <- NULL
    event_time <- NULL
    for(j in seq_along(Ynodes)) {
      if(data[i, Ynodes[j]] == 1) {
        first_event <- Ynodes[j]
        event_time <- ceiling(j/2)  # Convert to time period
        break
      }
    }
    
    if(!is.null(first_event)) {
      cat("\nExample row", i, "with first event at", first_event, ":\n")
      cat("Y nodes:", paste(Ynodes, "=", data[i, Ynodes], collapse=", "), "\n")
      cat("A nodes:", paste(Anodes, "=", data[i, Anodes], collapse=", "), "\n")
      cat("Treatment nodes after event should be NA\n")
    }
  }
}

# Run the check
check_data_patterns(data)

# Define node types
Anodes <- paste0("A", 1:3)
Ynodes <- c("D1", "Y1", "D2", "Y2", "D3", "Y3")
Lnodes <- NULL

#  deterministic Q function example
deterministic.Q.function <- function(data, current.node, nodes, called.from.estimate.g) {
  # Find death nodes by looking for 'D' followed by a number
  death_cols <- grep("^D[0-9]", names(data), value=TRUE)
  if(length(death_cols)==0) stop("No death/terminal event node found")
  
  # Get the Ynodes vector from the function environment
  Ynodes <- c("D1", "Y1", "D2", "Y2", "D3", "Y3")

  # Find position in node sequence
  current_name <- names(data)[current.node]
  if(is.null(current_name)) return(NULL)
  
  # If not a Y node, return NULL
  if(!(current_name %in% Ynodes)) return(NULL)
  
  # Find position in Ynodes
  current_node_position <- which(Ynodes == current_name)
  
  # Find death nodes that come before current node
  hist.death.index <- which(names(data) %in% death_cols[1:((current_node_position-1)/2)])
  
  if(length(hist.death.index)==0)
    return(NULL)
  else{
    is.deterministic <- Reduce("+",lapply(data[,hist.death.index,drop=FALSE],
                                          function(dd){x=dd;x[is.na(dd)] <- 0;x}))>=1
    is.deterministic[is.na(is.deterministic)] <- FALSE
    list(is.deterministic=is.deterministic, Q.value=1)
  }
}

# Run LTMLE
result_no_detQ <- ltmle(data,
                        Anodes = Anodes,
                        Cnodes = NULL,
                        Lnodes = Lnodes,
                        Ynodes = Ynodes,
                        survivalOutcome = TRUE,
                        abar = list(c(1,1,1), c(0,0,0)),
                deterministic.Q.function = NULL)
# Summarize results
summary(result_no_detQ)

# Run LTMLE without deterministic Q function
result_noDetQ <- ltmle(data,
                Anodes = Anodes,
                Cnodes = NULL,
                Lnodes = Lnodes,
                Ynodes = Ynodes,
                survivalOutcome = TRUE,
                abar = list(c(1,1,1), c(0,0,0)),
                deterministic.Q.function = NULL)


# Run LTMLE with deterministic Q function
result <- ltmle(data,
                Anodes = Anodes,
                Cnodes = NULL,
                Lnodes = Lnodes,
                Ynodes = Ynodes,
                survivalOutcome = TRUE,
                abar = list(c(1,1,1), c(0,0,0)),
                deterministic.Q.function = deterministic.Q.function)



# Look at patterns in the processed data
print("Sample of processed data:")
head(data[, Ynodes])

# Count events at each time
for(t in 1:3) {
  cat("\nTime", t, ":\n")
  cat("Deaths:", sum(data[[paste0("D", t)]]), "\n")
  cat("Outcomes:", sum(data[[paste0("Y", t)]]), "\n")
}


# Summarize results
summary(result_noDetQ)
summary(result)