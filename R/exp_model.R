exp_model <- function (family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea, marker, figure)
{

  #
  ssEnv <- get_session_info()
  exp_params <- unlist(strsplit(as.character(family_test),"_"))
  # log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " exp model fitting executing: ", as.character(sig.formula))
  partition_percentage <- as.numeric(exp_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)

  dep_var <- sig.formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

  if(length(exp_params)==3)
    if(exp_params[3]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }


  bias_dependent_variable <- min(tempDataFrame[, dependent_variable])
  tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1 - ifelse(bias_dependent_variable < 0, bias_dependent_variable, 0)
  # tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1

  # bias_independent_variable <- min(tempDataFrame[, independent_variable])
  # tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1 - ifelse(bias_independent_variable < 0, bias_independent_variable, 0)
  tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # Split the data into training and test set
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p=partition_percentage, list=FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # Assuming you have vectors train.data[,independent_variable] and train.data[,dependent_variable]
  start_a <- exp(stats::coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[1])  # Starting value for a, based on linear fit
  start_b <- coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[2]       # Starting value for b, based on linear fit

  exp_model_result <- NA
  tryCatch({
    stats::nls.control(warnOnly = TRUE)
    stats::nls.control(maxiter = 1000)
    stats::nls.control(minFactor = 1e-8)
    # Fit the exponential model
    formula <- as.formula(paste0(dependent_variable," ~ a * exp(b * ",independent_variable,")"))
    exp_model_result <- stats::nls(formula = formula , start = list(a = start_a, b = start_b), data = train.data)
  }, error = function(e) {

  })

  # check if exp_model_result is a model
  if (!inherits(exp_model_result, "nls"))
    return(res)

  predictions <- stats::predict(exp_model_result)
  if(nrow(test.data) != 0)
  {
    # Make predictions
    predictions_test <- stats::predict(exp_model_result, newdata = test.data)
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], c(),c()))

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(exp_model_result))
  # conf_int <- confint(exp_model_result)

  significative <- TRUE
  # for a and b extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- paste0("EXP_",row_name,"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name
    significative <- significative & p_value < as.numeric(ssEnv$alpha)

    exp_coef_estimate <- data.frame(exp_a_estimate = coefficients[i,1])
    colnames(exp_coef_estimate) <- toupper(paste0("EXP_",row_name,"_ESTIMATE",sep=""))

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

  if(plot & !any(is.na(predictions)))
  {
    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", name_cleaning(samples_sql_condition)))
    filename  =  file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),ssEnv$plot_format)
    #
    # ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
    #   ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
    #   ggplot2::stat_function( na.rm = TRUE, fun = function(x) (res$EXP_A_ESTIMATE * exp(res$EXP_B_ESTIMATE * x))  , color = ssEnv$color_palette_darker[3]) +
    #   ggplot2::ggtitle("") +
    #   ggplot2::xlab(independent_variable) +
    #   ggplot2::ylab(dependent_variable)

    ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_point(color = ssEnv$color_palette[1]) +
      ggplot2::stat_function(na.rm = TRUE, fun = function(x) (res$EXP_A_ESTIMATE * exp(res$EXP_B_ESTIMATE * x)), color = ssEnv$color_palette_darker[3]) +
      ggplot2::stat_smooth(method = "loess", color = ssEnv$color_palette[2], se = FALSE) +
      ggplot2::ggtitle("") +
      ggplot2::xlab(independent_variable) +
      ggplot2::ylab(dependent_variable)

    if(nrow(test.data) != 0)
    {
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3])+
        ggplot2::geom_point(data = test.data, ggplot2::aes(y = dependent_variable, x = eval(parse(text=independent_variable))), color = ssEnv$ssEnv$color_palette[2])+
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)
    }

    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 9,
      height = 9,
      units = c("in"),
      dpi = as.numeric(ssEnv$plot_resolution_ppi)
    )

    # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (start_a * exp(start_b * x))))
    # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")


  }
  return (res)
}

# exp_model_result("exp_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
