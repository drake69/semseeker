#' Title
#'
#' @param family_test family lqmm, quantreg
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param permutation_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
quantreg_model <- function(family_test, sig.formula, tempDataFrame, independent_variable, transformation, plot)
{

  # browser()
  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  res <- data.frame()
  tau <- as.numeric(quantreg_params[2])
  res <- data.frame("tau" = tau)

  dep_var <- sig.formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

  try({
    model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
    summary_qr <- suppressMessages(summary(model)$tTable)
    res$pvalue <- summary_qr[2,"Pr(>|t|)"]
    res$std.error <- summary_qr[2,"Std. Error"]
    res$regressor_coefficient <- summary_qr[2,"Value"]
    res$ci.lower <- summary_qr[2,"lower bound"]
    res$ci.upper <- summary_qr[2,"upper bound"]

    # predict.lqm
    # residuals.lqm
    predicted_values <- lqmm::predict.lqm(model, tempDataFrame)
    expected_values <- tempDataFrame[,dependent_variable]
    res <- quantreg_metrics(predicted_values = predicted_values, expected_values = expected_values,
      tau = tau, res = res, family_test = family_test, independent_variable = independent_variable, transformation = transformation,
      dependent_variable = dependent_variable, permutation_vector = c(), plot = plot)
  })
  res$r_model <- family_test

  #
  res <- as.data.frame(res)
  return (res)

}

