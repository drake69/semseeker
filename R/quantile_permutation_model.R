#' Title
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param tau tau at which apply the quantile regression
#' @param lqm_control specification of the lqmm package
#'
compute_quantile_delta_permutation <- function(sig.formula,df, shuffle = FALSE, quantile = 0.5)
{
  # #
  cols <- colnames(df)
  indepVar <- as.character(all.vars(sig.formula)[2])
  burden <- as.character(all.vars(sig.formula)[1])
  tempDataFrame <- df
  idx = which(colnames(tempDataFrame) == indepVar)
  tempDataFrame[,idx] <- as.factor(tempDataFrame[,idx])
  if (shuffle == TRUE)
    tempDataFrame[ ,indepVar] <- sample(tempDataFrame[,indepVar])
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  # calculate quantile difference based on indepVar
  quantile_difference <- tapply(tempDataFrame[,burden], tempDataFrame[,indepVar], quantile(probs=quantile))
  statistic_parameter <- diff(quantile_difference)
  statistic_parameter <- as.numeric(statistic_parameter[[1]])
  return(statistic_parameter)
}

#' quantile_permutation_model calculate differences between the same quantile of two distribution
#'
#' @param family_test family quantile
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param permutation_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
quantile_permutation_model <- function(family_test, sig.formula, tempDataFrame, independent_variable)
{
  quantile_params <- unlist(strsplit(as.character(family_test),"_"))

  # quantile_params template quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  # apply permutation to obtain signal using quantile
  # Define function to compute delta quantile regression coefficient

  n_permutations_test <- as.numeric(quantile_params[2])
  res <- data.frame(n_permutations = n_permutations_test)

  n_permutations <- as.numeric(quantile_params[3])
  conf.level <- as.numeric(quantile_params[4])
  res$conf.level <- conf.level
  quantile <- as.numeric(quantile_params[5])
  res$quantile <- quantile

  res$independent_variable <- as.character(all.vars(sig.formula)[2])

  pvalue_limit <- 1 - conf.level
  pvalue_limit_inf <- (pvalue_limit/2)
  pvalue_limit_sup <- 1 - (pvalue_limit/2)

  # Compute signal and p-value for n_permutations replications
  statistic_parameter <-  compute_quantile_delta_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle = FALSE)
  permutation_vector <- replicate(n_permutations_test, compute_quantile_delta_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle = TRUE, quantile = quantile))
  res$pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  res$statistic_parameter <- statistic_parameter
  if (pvalue>1)
    res$pvalue <- 1
  # Compute average signal and p-value
  if ((pvalue < pvalue_limit) && (n_permutations_test < n_permutations))
    permutation_vector <- replicate(n_permutations, compute_quantile_delta_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE, quantile = quantile))
  res$r_model <- "quantile_permutation_model"
  if(length(permutation_vector) == n_permutations_test)
    res$n_permutations <- n_permutations_test
  res$pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  ci <- stats::quantile(permutation_vector, probs = c(pvalue_limit_inf, pvalue_limit_sup))
  res$ci.lower <- ci[1]
  res$ci.upper <- ci[2]

  return (data.frame(ci.lower,ci.upper, pvalue, statistic_parameter,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations))
}
