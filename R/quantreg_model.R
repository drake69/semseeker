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

    # predicted_values <-predict(model, tempDataFrame)
    # res <- quantile_regression_metrics(predicted_values = model$fitted.values, expected_values = tempDataFrame[[independent_variable]],
    #   tau = tau, res = res, family_test = family_test, independent_variable = independent_variable, transformation = "quantile",
    #   dependent_variable = independent_variable, permutation_vector = c(), plot = FALSE)
  })
  res$r_model <- "lqmm_lqm"

  return (res)

}

# extract quantgile regression metrics and plot the results of the quantile regression and the permutation vector histogram
quantile_regression_metrics <- function(predicted_values, expected_values, tau, res, family_test,independent_variable,
  transformation, dependent_variable, permutation_vector=c(), plot = FALSE ){

  ###########

  ssEnv <- semseeker:::get_session_info()

  pinball_loss <- function(expected_values, predicted_values, tau) {
    residuals <- expected_values - predicted_values
    loss <- sum(ifelse(residuals < 0, (1 - tau) * -residuals, tau * residuals))
    return(loss)
  }

  # Calculate the Pinball Loss
  res$pinball_loss <- pinball_loss(expected_values, predicted_values, tau)

  ###########

  # Calculate the proportion of observations below the predicted quantile
  res$below_quantile <- mean(expected_values <= predicted_values)

  res$qq_distance <-  abs(tau - below_quantile)

  ###########

  if(plot)
  {
    chartFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  semseeker:::file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),"png")

    # Save the plot
    grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
    plot(expected_values, predicted_values, main = "Observed vs. Predicted Quantiles",
      xlab = "Observed Quantiles", ylab = "Predicted Quantiles", pch = 19)
    abline(0, 1, col = ssEnv$color_palette_darker[1] )  # Adding a 45-degree line for reference
    grDevices::dev.off()

    # plot histograms of betas
    if (length(permutation_vector) > 0)
    {
      chartFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
      filename  =  semseeker:::file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable, "HISTOGRAM"),"png")

      grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
      hist(expected_values - predicted_values, main = "Histogram of Betas",
        xlab = "Betas", ylab = "Frequency", col = ssEnv$color_palette[1])
      grDevices::dev.off()
    }
  }

  ###########

  # Calculate correlation
  correlation <- cor(expected_values, predicted_values)
  res$qq_correlation <- correlation

  ###########

  return(res)
}
