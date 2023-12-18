#' Title
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param tau tau at which apply the quantile regression
#' @param lqm_control specification of the lqmm package
#'
compute_mean_delta_boot <- function(sig.formula,df, shuffle = FALSE)
{
  # browser()
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
  # calculate mean difference based on indepVar
  mean_difference <- tapply(tempDataFrame[,burden], tempDataFrame[,indepVar], mean)
  beta_value <- diff(mean_difference)
  beta_value <- as.numeric(beta_value[[1]])
  return(beta_value)
}



#' Title
#'
#' @param family_test family lqmm, mean
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param boot_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
mean_bootstrap <- function(family_test, sig.formula, tempDataFrame, independent_variable, boot_success, tests_count)
{
  n_permutations <- NA
  ci.lower.adjusted <- NA
  ci.upper.adjusted <- NA
  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  beta_value <- NA
  pvalue <- NA
  boot.bca <- NA

  mean_params <- unlist(strsplit(as.character(family_test),"_"))

  # mean_params template mean + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  # apply bootstrap to obtain beta using mean
  # Define function to compute delta mean regression coefficient
  n_permutations_test <- as.numeric(mean_params[2])
  n_permutations <- as.numeric(mean_params[3])
  conf.level <- as.numeric(mean_params[4])
  conf.level = 1 - (1 - boot_success/tests_count) * (1 - conf.level)

  # Compute beta and p-value for n_permutations replications
  estimate <-  compute_mean_delta_boot(sig.formula=sig.formula, df=tempDataFrame, shuffle = FALSE)
  boot_vector <- replicate(n_permutations_test, compute_mean_delta_boot(sig.formula=sig.formula, df=tempDataFrame, shuffle = TRUE))
  p.value <- mean(abs(boot_vector) >= abs(estimate))
  # Compute average beta and p-value
  if ((p.value < 0.05) && (n_permutations_test < n_permutations))
  {
    boot_vector <- replicate(n_permutations, compute_mean_delta_boot(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE))
  }
  r_model <- "mean_bootstrap"
  boot.bca.adjusted <- coxed::bca(boot_vector, conf.level = conf.level)
  p.value <- mean(abs(boot_vector) >= abs(estimate))
  ci.lower.adjusted <-  boot.bca.adjusted[1]
  ci.upper.adjusted <- boot.bca.adjusted[2]

  return (data.frame(ci.lower,ci.upper, pvalue, beta_value,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations,ci.lower.adjusted,ci.upper.adjusted))

}
