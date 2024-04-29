#' Title
#'
#' @param family_test family lqmm, quantreg
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param permutation_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
quantreg_model <- function(family_test, sig.formula, tempDataFrame, independent_variable)
{

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  res <- data.frame()
  tau = as.numeric(quantreg_params[2])
  res <- data.frame("tau" = tau)
  try({
    model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
    summary_qr <- suppressMessages(summary(model)$tTable)
    res$pvalue <- summary_qr[2,"Pr(>|t|)"]
    res$std.error <- summary_qr[2,"Std. Error"]
    res$statistic_parameter <- summary_qr[2,"Value"]
    res$ci.lower <- summary_qr[2,"lower bound"]
    res$ci.upper <- summary_qr[2,"upper bound"]
  })
  res$r_model <- "lqmm_lqm"

  return (res)

}
