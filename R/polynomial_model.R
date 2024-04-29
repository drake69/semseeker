polynomial_model <- function (family_test, tempDataFrame, sig.formula , transformation, plot)
{

  ssEnv <- semseeker:::get_session_info()

  # plynomial_degree_partition-partition_percentage
  polynomial_params <- unlist(strsplit(as.character(family_test),"_"))

  degree <- as.numeric(polynomial_params[2])
  res <- data.frame("PL_DEGREE"= degree)
  partition_percentage <- as.numeric(polynomial_params[3])
  res$PL_PERC <- partition_percentage
  # res$sig.fromula <- as.character(sig.formula)

  if(length(polynomial_params)==4)
    res$r_model <- paste0("polynomial",polynomial_params[4], sep="")
  else
    res$r_model <- "polynomial"


  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]
  if(length(polynomial_params)==4)
    if(polynomial_params[4]=="predictor")
    {
      dependent_variable <- dep_var[[3]]
      independent_variable <- dep_var[[2]]
    }

  # Split the data into training and test set
  training.samples <- tempDataFrame[, dependent_variable] %>% caret::createDataPartition(p = partition_percentage, list = FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # Build the model
  model <- stats::lm(eval(parse(text=dependent_variable)) ~ stats::poly(eval(parse(text=independent_variable)), degree, raw = TRUE), data = train.data)

  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(model$residuals^2)
  # Calculate Root Mean Square Error (RMSE)
  res$rmse <- sqrt(mean(model$residuals^2))
  # Goodness-of-Fit Statistics (R-squared)
  res$r_squared <- summary(model)$r.squared
  # Calculate Adjusted R-squared
  res$adj_r_squared <- summary(model)$adj.r.squared
  # Mean Absolute Error (MAE)
  res$mae <- mean(abs(model$residuals))
  # Mean Squared Log Error (MSLE)
  res$msle <- mean((log(predict(model) + 1) - log(train.data[, dependent_variable] + 1))^2)

  if(nrow(test.data) != 0)
  {
    # Make predictions
    predictions_test <- stats::predict(model, newdata = test.data)
    # Calculate RMSE test
    res$rmse_test = caret::RMSE(predictions_test, test.data[,dependent_variable])
    # Calculate Root Mean Square Error (RMSE)
    res$r_squared_test = caret::R2(predictions_test, test.data[,dependent_variable])
  }

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(model))
  # conf_int <- confint(model)

  significative <- TRUE
  # for each degree extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    if (row_name =="(Intercept)")
      pval_name <- "PL_INTERCEPT_PVALUE"
    else
      pval_name <- paste0("PL_DEGREE_",as.character(i -1 ),"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    significative <- significative & p_value < ssEnv$alpha
    colnames(p_value) <- pval_name
    res <- cbind(res, p_value)
  }

  res$pvalue <- max(coefficients[,4])
  colnames(significative) <- "SIGNIFICATIVE"
  res <- cbind(res, significative)
  # remove rowname from res
  rownames(res) <- NULL
  if(plot)
  {
    chartFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  semseeker:::file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),"png")
    # Plotting the fitted curve
    # ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
    #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
    #   ggplot2::stat_function(fun = function(x) start_a * exp(start_b * x), color = ssEnv$color_palette[1]) +
    #   ggplot2::ggtitle("Fitted Exponential Curve") +
    #   ggplot2::xlab(independent_variable) +
    #   ggplot2::ylab(dependent_variable)

    # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
    ggp <- ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
      ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
      ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE), color = ssEnv$color_palette_darker[2])

    if (nrow(test.data) != 0)
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3]) +
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
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
  #   ggplot2::stat_smooth(method = lm, formula = y ~ stats::poly(x, degree, raw = TRUE))
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2])
  #
  # # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
  # ggplot2::ggplot(train.data, ggplot2::aes(eval(parse(text=independent_variable)), eval(parse(text=dependent_variable))) ) +
  #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) + ggplot2::stat_smooth(method = lm, formula = y ~ poly(x, degree, raw = TRUE)) +
  #   ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions), color = ssEnv$color_palette[2]) +
  #   ggplot2::geom_point(data = data.frame(train.data,model$residuals) , ggplot2::aes(y = model$residuals), color = "cyan")
  #
  return (res)
}

# polynomial_model("polynomial_4_0.8", sample_sheet, "DELTARQ_HYPO ~ Age")
