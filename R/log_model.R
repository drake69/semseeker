log_model <- function (family_test, tempDataFrame, sig.formula)
{
  # plynomial_degree_partition-partition_percentage
  log_params <- unlist(strsplit(as.character(family_test),"_"))

  partition_percentage <- as.numeric(log_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  if(length(log_params)==3)
    res$r_model <- paste0("log",log_params[3], sep="")
  else
    res$r_model <- "log"


  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]
  if(length(log_params)==4)
    if(log_params[4]=="predictor")
    {
      dependent_variable <- dep_var[[3]]
      independent_variable <- dep_var[[2]]
    }

  tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1
  tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # Split the data into training and test set
  training.samples <- tempDataFrame[, dependent_variable] %>% caret::createDataPartition(p = partition_percentage, list = FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]


  x_train_data <- train.data[, independent_variable]
  y_train_data <- train.data[, dependent_variable]
  # Fit the logarithmic log_model
  log_model <- stats::lm((y_train_data) ~ log(x_train_data))

  # Calculate Adjusted R-squared
  res$adj_r_squared <- summary(log_model)$adj.r.squared
  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(log_model$residuals^2)
  # Calculate Root Mean Square Error (RMSE)
  res$rmse <- sqrt(mean((predict(log_model) - (y_train_data))^2))
  # Goodness-of-Fit Statistics (R-squared)
  res$r_squared <- summary(log_model)$r.squared


  # Make predictions
  predictions <- stats::predict(log_model)
  # log_model performance
  res$rmse_test = caret::RMSE(predictions, test.data[,dependent_variable])
  # res$r_squared_test = caret::R2(predictions, test.data[,dependent_variable])

  # browser()
  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(log_model))
  # conf_int <- confint(log_model)

  # for a and b extract the p-value
  # for (i in 1:nrow(coefficients)) {
  #   # i <- 1
  #   p_value <- coefficients[i,4]
  #   row_name <- rownames(coefficients)[i]
  #   pval_name <- paste0("log_",row_name,"_PVALUE",sep="")
  #   p_value <- data.frame(p_value)
  #   colnames(p_value) <- pval_name
  #
  #   log_coef_estimate <- data.frame(log_a_estimate = coefficients[i,1])
  #   colnames(log_coef_estimate) <- paste0("log_",row_name,"_ESTIMATE",sep="")
  #
  #   res <- cbind(res, p_value, log_coef_estimate)
  # }

  res$pvalue <- coefficients[2,4]
  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  return (res)
}

# log_model("log_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
