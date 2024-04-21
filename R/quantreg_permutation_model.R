#' Title
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param tau tau at which apply the wuantile regression
#' @param lqm_control specification of the lqmm package
#'
compute_quantreg_permutation <- function(sig.formula,df, tau, lqm_control)
{
  # #
  cols <- colnames(df)
  # tempDataFrame <- Rfast::colShuffle(as.matrix(df))
  indepVar <- as.character(all.vars(sig.formula)[-1][1])
  tempDataFrame <- df
  tempDataFrame[ ,indepVar] <- sample(tempDataFrame[,indepVar])
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  suppressMessages({
    model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
  })
  summary_qr <- suppressMessages(summary(model)$tTable)
  not_permutated_regression_coefficient <- summary_qr[2,"Value"]
  return(not_permutated_regression_coefficient)
}

#' Title
#'
#' @param family_test family lqmm, quantreg
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param permutation_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
quantreg_permutation_model <- function(family_test, sig.formula, tempDataFrame, independent_variable)
{
  n_permutations <- NA
  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  not_permutated_regression_coefficient <- NA
  pvalue <- NA
  summary_results <- NA

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))

  # do permutations
  if(length(quantreg_params)!=5)
  {
    log_event("ERROR: number of parameter incorrect for quantreg-permutation expected:
      quantreg-permutation + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_regression_coefficient")
    exit()
  }

  tau <- as.numeric(quantreg_params[2])
  n_permutations_test <- as.numeric(quantreg_params[3])
  n_permutations <- as.numeric(quantreg_params[4])
  conf.level <- as.numeric(quantreg_params[5])

  pvalue_limit <- 1 - conf.level
  pvalue_limit_inf <- (pvalue_limit/2)
  pvalue_limit_sup <- 1 - (pvalue_limit/2)

  suppressMessages({
    model <- lqmm::lqm(sig.formula, tau =tau, data = as.data.frame(tempDataFrame), na.action = stats::na.omit, control = lqm_control)
  })
  summary_qr <- suppressMessages(summary(model)$tTable)
  not_permutated_regression_coefficient <- summary_qr[2,"Value"]
  std.error <- summary_qr[2,"Std. Error"]
  ci.lower <- summary_qr[2,"lower bound"]
  ci.upper <- summary_qr[2,"upper bound"]
  pvalue <- summary_qr[2,"Pr(>|t|)"]

  if(n_permutations > n_permutations_test)
  {
    permutated_regression_coefficient <- suppressMessages(replicate(n_permutations_test, compute_quantreg_permutation(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control)))
    summary_results <- exact_pvalue(permutated_regression_coefficient, not_permutated_regression_coefficient, conf.level = conf.level)
    pvalue <- summary_results[3]
  }
  if(pvalue < pvalue_limit)
  {
    permutated_regression_coefficient <- suppressMessages(replicate(n_permutations, compute_quantreg_permutation(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control)))
    summary_results <- exact_pvalue(permutated_regression_coefficient, not_permutated_regression_coefficient, conf.level = conf.level)
    pvalue <- summary_results[3]
  }
  ci.lower <- summary_results[1]
  ci.upper <- summary_results[2]

  if(length(permutated_regression_coefficient) == n_permutations_test)
    n_permutations <- n_permutations_test

  r_model <- "lqmm_lqm"

  return (data.frame(ci.lower,ci.upper, pvalue, not_permutated_regression_coefficient,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations))

}
