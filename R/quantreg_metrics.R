# extract quantgile regression metrics and plot the results of the quantile regression and the permutation vector histogram
quantreg_metrics <- function(predicted_values, expected_values, tau, res, family_test,independent_variable,
  transformation, dependent_variable, permutation_vector=c(), plot = FALSE ){

  ###########

  ssEnv <- get_session_info()

  pinball_loss <- function(expected_values, predicted_values, tau) {
    residuals <- expected_values - predicted_values
    loss <- sum(ifelse(residuals < 0, (1 - tau) * (-residuals), tau * residuals))
    return(loss)
  }

  # Calculate the Pinball Loss
  local_res <- data.frame("pinball_loss" = pinball_loss(expected_values, predicted_values, tau))

  ###########

  # Calculate the proportion of observations below the predicted quantile
  local_res$below_quantile <- mean(expected_values <= predicted_values)

  local_res$qq_distance <-  abs(tau - local_res$below_quantile)

  ###########

  if(plot)
  {
    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),ssEnv$plot_format)

    # Save the plot
    # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
    if(ssEnv$plot_format == "png")
      grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300, bg = "transparent")
    if(ssEnv$plot_format == "eps")
      grDevices::postscript(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300, bg = "transparent")
    graphics::plot(expected_values, predicted_values, main = "Observed vs. Predicted",
      xlab = "Observed", ylab = "Predicted", pch = 19)
    graphics::abline(0, 1, col = ssEnv$color_palette_darker[1] )  # Adding a 45-degree line for reference
    grDevices::dev.off()

    # plot histograms of betas
    if (length(permutation_vector) > 0)
    {
      chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
      filename  =  file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable, "HISTOGRAM"),ssEnv$plot_format)

      # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
      if(ssEnv$plot_format == "png")
        grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300, bg = "transparent")
      if(ssEnv$plot_format == "eps")
        grDevices::postscript(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300, bg = "transparent")
      graphics::hist(permutation_vector, main = "Histogram of Betas",
        xlab = "Regression Coefficent Value", ylab = "Frequency", col = ssEnv$color_palette[1])
      grDevices::dev.off()
    }
  }

  ###########

  # Calculate correlation
  correlation <- stats::cor(expected_values, predicted_values)
  local_res$qq_correlation <- correlation

  ###########
  model_performance_res <- model_performance(predicted_values, expected_values, c(),c())
  res <- cbind(res, local_res, model_performance_res)

  return(res)
}
