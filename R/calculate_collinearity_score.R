# Function to calculate the collinearity score among covariates
calculate_collinearity_score <- function(df) {
  # # Install and load required packages
  # if (!require("car")) install.packages("car")
  # library(car)
  # Calculate VIF for each variable using a linear model where each variable is regressed on all others

  # browser()
  vif_values <- apply(df, 2, function(var) {
    lm_model <- stats::lm(var ~ ., data = df)
    return(max(car::vif(lm_model)))
  })

  # Calculate Condition Index
  X <- as.matrix(df)
  condition_index <- kappa(X, exact = TRUE)

  # Perform PCA
  pca <- stats::prcomp(df, scale. = TRUE)
  explained_variance <- summary(pca)$importance[2, ]

  # Normalizing and weighting scores
  vif_score <- mean(vif_values)
  condition_index_score <- condition_index / 30  # Assuming 30 as a threshold
  pca_score <- sum(explained_variance[explained_variance > 0.1])  # Components explaining >10% variance

  # Aggregate score (weights can be adjusted as needed)
  total_score <- (vif_score * 0.5) + (condition_index_score * 0.3) + (pca_score * 0.2)

  if((vif_score > 5) | (condition_index_score > 0.5) | (pca_score < 0.5)) {
    log_event("ERROR: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Collinearity check failed. VIF: ", vif_score, ", Condition Index: ", condition_index_score, ", PCA: ", pca_score)
  }

  return(total_score)
}

