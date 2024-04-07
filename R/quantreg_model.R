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
  n_permutations <- NA
  
  
  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  statistic_parameter <- NA
  pvalue <- NA
  

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  tau = as.numeric(quantreg_params[2])
  try({
    model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
    summary_qr <- suppressMessages(summary(model)$tTable)
    pvalue <- summary_qr[2,"Pr(>|t|)"]
    std.error <- summary_qr[2,"Std. Error"]
    statistic_parameter <- summary_qr[2,"Value"]
    ci.lower <- summary_qr[2,"lower bound"]
    ci.upper <- summary_qr[2,"upper bound"]
  })
  r_model <- "lqmm_lqm"

  return (data.frame(ci.lower,ci.upper, pvalue, statistic_parameter,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations))

}
