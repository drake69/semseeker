#
# # Install and load necessary packages
# if (!require("glmnet::glmnet")) install.packages("glmnet::glmnet")
# if (!require("boot")) install.packages("boot")
#
# library(glmnet::glmnet)
# library(boot)

# Function to calculate mediation effects using Ridge Regression
mediation_ridge <-  function(family_test,tempDataFrame, sig.formula, transformation, plot) {
  # function(tempDataFrame, treatment, mediator, outcome,covariates, lambda = 0.1, permutations = 1000) {

  # Ensure columns exist
  if (!all(c(treatment, mediator, outcome) %in% colnames(tempDataFrame))) {
    stop("One or more specified columns do not exist in the tempDataFrame.")
  }

  ssEnv <- get_session_info()

  ridge_params <- unlist(strsplit(as.character(family_test),"_"))
  lambda <- as.numeric(ridge_params[2])
  res <- data.frame("lambda" = lambda)
  permutations_test <- as.numeric(ridge_params[3])
  permutations <- as.numeric(ridge_params[4])

  # decompose the formula, to move the mediator
  # the expected formula is "mediator ~ outcome + treatment + covariates"
  vars <- sig.formula_vars(sig.formula)
  mediator <- vars[[1]] # burden_Value
  outcome <- vars[[2]] # dependent variable
  treatment <- vars[[3]][1]

  log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Starting mediation analysis for: ", mediator, " and ", outcome)

  res$mediator <- mediator
  res$outcome <- outcome
  res$treatment <- treatment

  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- sapply(covariates, as.character)
  } else
    covariates <- c()


  # Prepare data matrices for glmnet::glmnet
  if(length(covariates) != 0)
    x <- as.matrix(tempDataFrame[, c(treatment, covariates)])
  else
    x <- as.matrix(tempDataFrame[, c(treatment)])
  z <- as.matrix(tempDataFrame[, mediator])
  y <- as.matrix(tempDataFrame[, outcome])

  # Ensure the matrices are correctly dimensioned
  x <- matrix(x, ncol = (1 + length(covariates)))
  z <- matrix(z, ncol = 1)
  y <- matrix(y, ncol = 1)

  perms <- sort(unique(c(permutations_test, permutations)))

  # Fit the mediator model using Ridge Regression
  mediator_model <- glmnet::glmnet(x, z, alpha = 0, lambda = lambda)
  treatment_on_mediator_effect <- coef(mediator_model)[2]

  # Fit the outcome model using Ridge Regression
  outcome_model <- glmnet::glmnet(cbind(x, z), y, alpha = 0, lambda = lambda)
  mediator_on_outcome_effect <- coef(outcome_model)[3]
  direct_effect <- coef(outcome_model)[2]

  # Calculate the indirect effect (indirect_effect)
  indirect_effect <- treatment_on_mediator_effect * mediator_on_outcome_effect

  # Calculate the total effect (c = c' + indirect_effect)
  total_effect <- direct_effect + indirect_effect

  browser()
  for (i in 1:length(perms)) {

    permutations <- perms[i]

    # Perform bootstrapping
    boot_results <- boot::boot(data = tempDataFrame, statistic = boot_mediation_ridge, R = permutations)

    # Calculate p-values
    p_value <- function(observed, bootstrapped) {
      (sum(abs(bootstrapped) >= abs(observed)) + 1) / (length(bootstrapped) + 1)
    }

    p_values <- sapply(1:4, function(i) p_value(c(treatment_on_mediator_effect, mediator_on_outcome_effect, direct_effect, indirect_effect)[i], boot_results$t[, i]))
    p_values <- as.data.frame(matrix(p_values, nrow = 1))
    colnames(p_values) <- c("Treatment_on_Mediator_Pvalue", "Mediator_on_Outcome_Pvalue", "Direct_Effect_Pvalue", "Indirect_Effect_Pvalue")

    # Calculate confidence intervals
    mediation_res <- boot_results[1,]
    boot.res <- boot::boot.ci(mediation_res, type="bca")
    res$treatment_on_mediator_effect_ci_lower <- boot.res$bca[4]
    res$treatment_on_mediator_effect_ci_upper <- boot.res$bca[5]
    # res$treatment_on_mediator_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

    # get proportion of the effect
    prop_res <- boot_results[2,]
    boot.res <- boot::boot.ci(prop_res, type="bca")
    res$mediator_on_outcome_effect_ci_lower <- boot.res$bca[4]
    res$mediator_on_outcome_effect_ci_upper <- boot.res$bca[5]
    # res$mediator_on_outcome_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

    #get direct effect
    prop_res <- boot_results[3,]
    boot.res <- boot::boot.ci(prop_res, type="bca")
    res$direct_effect_ci_lower <- boot.res$bca[4]
    res$direct_effect_ci_upper <- boot.res$bca[5]
    # res$direct_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

    #get indirect effect
    prop_res <- boot_results[4,]
    boot.res <- boot::boot.ci(prop_res, type="bca")
    res$indirect_effect_ci_lower <- boot.res$bca[4]
    res$indirect_effect_ci_upper <- boot.res$bca[5]
    # res$indirect_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

    #get total effect
    prop_res <- boot_results[5,]
    boot.res <- boot::boot.ci(prop_res, type="bca")
    res$total_effect_ci_lower <- boot.res$bca[4]
    res$total_effect_ci_upper <- boot.res$bca[5]
    # res$total_effect_significative <- !(boot.res$bca[4]<0 & boot.res$bca[5]>0)

    # # Create treatment_on_mediator_effect list to return the results
    # res <- list(
    #   treatment_on_mediator_effect = treatment_on_mediator_effect,
    #   mediator_on_outcome_effect = mediator_on_outcome_effect,
    #   direct_effect = direct_effect,
    #   indirect_effect = indirect_effect,
    #   total_effect = total_effect,
    #   p_values = p_values,
    #   boot_results = boot_results
    # )

    res <- cbind(res, pvalue)
    if(!all(pvalues < as.numeric(ssEnv$alpha)))
      break

  }
  return(res)
}

# # Example usage
# set.seed(123)
# df <- data.frame(
#   treatment = rbinom(100, 1, 0.5),
#   mediator = rnorm(100),
#   outcome = rnorm(100),
#   covariates = rnorm(100)
# )
#
# # Calculate mediation effects with Ridge Regression and bootstrapping
# mediation_result <- calculate_mediation_effect_ridge(df, treatment = "treatment", mediator = "mediator", outcome = "outcome", covariates = "covariates")
#
# # Print the results
# print(mediation_result)

# Function to calculate mediation effects for bootstrapped samples
boot_mediation_ridge <- function(tempDataFrame, indices) {
  boot_data <- tempDataFrame[indices, ]
  if(length(covariates) != 0)
    x_boot <- as.matrix(boot_data[, c(treatment, covariates)])
  else
    x_boot <- as.matrix(boot_data[, c(treatment)])
  z_boot <- as.matrix(boot_data[, mediator])
  y_boot <- as.matrix(boot_data[, outcome])

  x_boot <- matrix(x_boot, ncol = (1 + length(covariates)))
  z_boot <- matrix(z_boot, ncol = 1)
  y_boot <- matrix(y_boot, ncol = 1)

  mediator_model <- glmnet::glmnet(x_boot, z_boot, alpha = 0, lambda = lambda)
  outcome_model <- glmnet::glmnet(cbind(x_boot, z_boot), y_boot, alpha = 0, lambda = lambda)
  treatment_on_mediator_effect <- coef(mediator_model)[2]
  mediator_on_outcome_effect <- coef(outcome_model)[3]
  direct_effect <- coef(outcome_model)[2]
  indirect_effect <- treatment_on_mediator_effect * mediator_on_outcome_effect
  total_effect <- direct_effect + indirect_effect
  return(c(treatment_on_mediator_effect, mediator_on_outcome_effect, direct_effect, indirect_effect, total_effect))
}
