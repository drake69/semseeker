log_model <- function (family_test, tempDataFrame, sig.formula, transformation, plot)
{

  ssEnv <- semseeker:::get_session_info()
  # plynomial_degree_partition-partition_percentage
  log_params <- unlist(strsplit(as.character(family_test),"_"))
  partition_percentage <- as.numeric(log_params[2])
  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]

  tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1
  tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # Split the data into training and test set
  training.samples <- tempDataFrame[, dependent_variable] %>% caret::createDataPartition(p = partition_percentage, list = FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  train.data[, dependent_variable] <- train.data[, dependent_variable] + 1
  formula <- as.formula(paste0(dependent_variable, "~ log(", independent_variable, ")"))
  # Fit the logarithmic log_model_result
  log_model_result <- stats::lm(formula = formula, data = train.data)

  # if log_model_result a model
  if (is.null(log_model_result))
    return(res)

  predictions <- stats::predict(log_model_result)
  if(nrow(test.data) != 0)
  {
    # Make predictions
    predictions_test <- stats::predict(log_model_result, newdata = test.data)
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], c(),c()))

  if(plot)
  {
    chartFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  semseeker:::file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),"png")

    # do a plot with train.data, test.data and predictions with 3 different colors 1 color for train.data, 1 color for test.data and 1 color for predictions
    ggp <- ggplot2::ggplot(train.data, ggplot2::aes(((independent_variable)), ((dependent_variable))) ) +
      ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
      ggplot2::stat_smooth(method = lm, formula = y ~ log(x), se = FALSE, color = ssEnv$color_palette_darker[2])

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

  # browser()
  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(log_model_result))
  # conf_int <- confint(log_model_result)

  res$pvalue <- coefficients[2,4]
  # remove rowname from res
  rownames(res) <- NULL

  # plot log fitted log_model_result



  return (res)
}

# log_model_result("log_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
