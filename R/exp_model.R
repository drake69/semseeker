exp_model <- function (family_test, tempDataFrame, sig.formula, transformation, plot)
{

  ssEnv <- get_session_info()
  exp_params <- unlist(strsplit(as.character(family_test),"_"))

  partition_percentage <- as.numeric(exp_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]

  bias <- min(tempDataFrame[, dependent_variable])
  if (bias < 0)
    tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] - bias +1
  else
    tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1

  tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # Split the data into training and test set
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p=partition_percentage, list=FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # Assuming you have vectors train.data[,independent_variable] and train.data[,dependent_variable]
  start_a <- exp(stats::coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[1])  # Starting value for a, based on linear fit
  start_b <- coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[2]       # Starting value for b, based on linear fit

  # Fit the exponential model
  formula <- as.formula(paste0(dependent_variable," ~ a * exp(b * ",independent_variable,")"))
  exp_model <- nls( formula = formula , start = list(a = start_a, b = start_b), data = train.data)

  predictions <- stats::predict(exp_model)
  # Calculate residuals
  model_residuals <- predictions - train.data[, dependent_variable]
  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(model_residuals^2)
  # Calculate Root Mean Square Error (RMSE)
  res$rmse <- sqrt(mean((model_residuals)^2))
  # Goodness-of-Fit Statistics (R-squared)
  res$r_squared = caret::R2(predictions, train.data[,dependent_variable])
  # Calculate Mean Absolute Error (MAE)
  res$mae <- mean(abs(model_residuals))
  # Calculate Mean Squared Log Error (MSE)
  res$msle <- mean((log(predictions + 1) - log(train.data[, dependent_variable] + 1))^2)

  if(nrow(test.data) != 0)
  {
    # Make predictions
    predictions_test <- stats::predict(exp_model, newdata = test.data)
    # Calculate RMSE test
    res$rmse_test = caret::RMSE(predictions_test, test.data[,dependent_variable])
    # Calculate Root Mean Square Error (RMSE)
    res$r_squared_test = caret::R2(predictions_test, test.data[,dependent_variable])
  }

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(exp_model))
  # conf_int <- confint(exp_model)

  significative <- TRUE
  # for a and b extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- paste0("EXP_",row_name,"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name
    significative <- significative & p_value < ssEnv$alpha

    exp_coef_estimate <- data.frame(exp_a_estimate = coefficients[i,1])
    colnames(exp_coef_estimate) <- paste0("EXP_",row_name,"_ESTIMATE",sep="")

    res <- cbind(res, p_value, exp_coef_estimate)
  }

  # Add the max p-value
  res$pvalue <- max(coefficients[,4])
  # Add the significative
  colnames(significative) <- "SIGNIFICATIVE"
  res <- cbind(res, significative)
  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  if(plot)
  {

    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),"png")
    # Plotting the fitted curve
    ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
      ggplot2::stat_function(fun = function(x) start_a * exp(start_b * x), color = ssEnv$color_palette_darker[2]) +
      ggplot2::ggtitle("") +
      ggplot2::xlab(independent_variable) +
      ggplot2::ylab(dependent_variable)

    if(nrow(test.data) != 0)
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3])+
            ggplot2::xlab(independent_variable) +
            ggplot2::ylab(dependent_variable)


    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 1240,
      height = 1240,
      units = c("px"),
      dpi = 144
    )
  }
  return (res)
}

# exp_model("exp_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
