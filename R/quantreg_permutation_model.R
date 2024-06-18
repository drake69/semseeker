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
quantreg_permutation_model <- function(family_test, sig.formula, tempDataFrame, independent_variable, transformation, plot)
{
  #
  ssEnv <- get_session_info()

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))

  # do permutations
  if(length(quantreg_params)!=5)
  {
    log_event("ERROR: number of parameter incorrect for quantreg-permutation expected:
      quantreg-permutation + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_regression_coefficient")
    stop()
  }

  dep_var <- sig.formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

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
  res <- data.frame("not_permutated_regression_coefficient" = not_permutated_regression_coefficient)
  res$std.error <- summary_qr[2,"Std. Error"]
  res$ci.lower <- summary_qr[2,"lower bound"]
  res$ci.upper <- summary_qr[2,"upper bound"]
  res$pvalue <- summary_qr[2,"Pr(>|t|)"]
  res$conf.level <- conf.level

  perms <- sort(unique(c(n_permutations, n_permutations_test)), decreasing = F)

  for(p in 1:length(perms))
  {
    perm <- perms[p]
    res$n_permutations <- perm
    permutated_regression_coefficient <- suppressMessages(replicate(perm, compute_quantreg_permutation(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control)))
    summary_results <- exact_pvalue(permutated_regression_coefficient, not_permutated_regression_coefficient, conf.level = conf.level)
    res$pvalue <- summary_results[3]
    if(res$pvalue > ssEnv$alpha)
      break
  }

  res$ci.lower <- summary_results[1]
  res$ci.upper <- summary_results[2]

  predicted_values <- lqmm::predict.lqm(model, tempDataFrame)
  expected_values <- tempDataFrame[,dependent_variable]
  res <- quantreg_metrics(predicted_values = predicted_values, expected_values = expected_values,
    tau = tau, res = res, family_test = family_test, independent_variable = independent_variable, transformation = transformation,
    dependent_variable = dependent_variable, permutation_vector = permutated_regression_coefficient, plot = plot)

  if(length(permutated_regression_coefficient) == n_permutations_test)
    n_permutations <- n_permutations_test

  res$n_permutations <- n_permutations
  res$r_model <- family_test

  return (res)

}
