model_performance <- function(fitted_values, expected_values, prediction_values, prediction_expected_values)
{

  model_residuals <- expected_values - fitted_values

  # calculate the mean squared error
  mse <- mean((model_residuals)^2)
  res <- data.frame("mse"=mse)

  # calculate the root mean squared error
  res$rmse <- sqrt(mse)

  # calculate the mean absolute percentage error
  res$mape <- mean(abs(model_residuals/fitted_values), na.rm = TRUE)

  # calculate the mean percentage error
  res$mpe <- mean(model_residuals/fitted_values, na.rm = TRUE)

  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(model_residuals^2)

  # Calculate Mean Absolute Error (MAE)
  res$mae <- mean(abs(model_residuals), na.rm = TRUE)

  # calculate r-squared the higher the better - goodness of fit
  res$r_squared <- 1 - res$sse / sum((expected_values - mean(expected_values))^2)
  # complement of r-squared the lower the better
  res$r_squared_compl <- res$sse / sum((expected_values - mean(expected_values))^2)

  # calculate adjusted r-squared
  res$r_squared_adj <- 1 - (1 - res$r_squared) * (length(expected_values) - 1) / (length(expected_values) - length(fitted_values) - 1)
  # complement of adjusted r-squared the lower the better
  res$r_squared_adj_compl <- (1 - res$r_squared) * (length(expected_values) - 1) / (length(expected_values) - length(fitted_values) - 1)

  # Calculate Mean Squared Log Error (MSLE) the lower the better
  res$msle <- mean((log(fitted_values + 1) - log(expected_values + 1))^2, na.rm = TRUE)


  if(length(prediction_values)>0)
  {
    model_residuals <- prediction_expected_values - prediction_values
    expected_values <- prediction_expected_values
    fitted_values <- prediction_values

    # calculate the mean squared error
    res$mse_test <- mean((model_residuals)^2, na.rm = TRUE)

    # calculate the root mean squared error
    res$rmse_test <- sqrt(res$mse_test)

    # calculate the mean absolute percentage error
    res$mape_test <- mean(abs(model_residuals/fitted_values), na.rm = TRUE)

    # calculate the mean percentage error
    res$mpe_test <- mean(model_residuals/fitted_values)

    # Calculate Sum of Squares due to Error (SSE)
    res$sse_test <- sum(model_residuals^2)

    # Calculate Mean Absolute Error (MAE)
    res$mae_test <- mean(abs(model_residuals), na.rm = TRUE)

    # calculate r-squared the higher the better - goodness of fit
    res$r_squared_test <- 1 - res$sse_test / sum((expected_values - mean(expected_values))^2)
    # complement of r-squared the lower the better
    res$r_squared_compl_test <- res$sse_test / sum((expected_values - mean(expected_values))^2)

    # calculate adjusted r-squared
    res$r_squared_adj_test <- 1 - (1 - res$r_squared_test) * (length(expected_values) - 1) / (length(expected_values) - length(fitted_values) - 1)
    # complement of adjusted r-squared the lower the better
    res$r_squared_adj_compl_test <- (1 - res$r_squared_test) * (length(expected_values) - 1) / (length(expected_values) - length(fitted_values) - 1)

    # Calculate Mean Squared Log Error (MSLE) the lower the better
    res$msle_test <- mean((log(fitted_values + 1) - log(expected_values + 1))^2, na.rm = TRUE)

    if (rse$mse > rse$mse_test | rse$rmse > rse$rmse_test |
        rse$mape > rse$mape_test | rse$mpe > rse$mpe_test |
        rse$mae > rse$mae_test | rse$r_squared < rse$r_squared_test |
        rse$r_squared_adj < rse$r_squared_adj_test |
        rse$msle > rse$msle_test)
      res$overfitting <- TRUE
    else
      res$overfitting <- FALSE
  }

  return(res)
}
