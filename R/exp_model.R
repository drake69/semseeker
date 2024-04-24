exp_model <- function (family_test, tempDataFrame, sig.formula)
{

  # plynomial_degree_partition-partition_percentage
  exp_params <- unlist(strsplit(as.character(family_test),"_"))

  partition_percentage <- as.numeric(exp_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  if(length(exp_params)==3)
    res$r_model <- paste0("exp",exp_params[3], sep="")
  else
    res$r_model <- "exp"


  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]
  if(length(exp_params)==4)
    if(exp_params[4]=="predictor")
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

  # Assuming you have vectors train.data[,independent_variable] and train.data[,dependent_variable]
  start_a <- exp(stats::coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[1])  # Starting value for a, based on linear fit
  start_b <- coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[2]       # Starting value for b, based on linear fit

  x_train_data <- train.data[, independent_variable]
  y_train_data <- train.data[, dependent_variable]
  # Fit the exponential model
  model <- nls(y_train_data ~ a * exp(b * x_train_data), start = list(a = start_a, b = start_b), data = data.frame(x_train_data, y_train_data))

  # Build the model
  # model <- stats::lm(eval(parse(text=dependent_variable)) ~ stats::poly(eval(parse(text=independent_variable)), degree, raw = TRUE), data = train.data)
  # model <- stats::lm(train.data[,dependent_variable] ~ stats::poly(train.data[,independent_variable], degree, raw = TRUE))
  # dependent_variable_train <- train.data[,dependent_variable]
  # independent_variable_train <- train.data[,independent_variable]
  # model <- stats::lm(dependent_variable_train ~ stats::poly(independent_variable_train, degree, raw = TRUE))

  # Calculate Adjusted R-squared
  res$adj_r_squared <- summary(model)$adj.r.squared
  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(model$residuals^2)
  # Calculate Root Mean Square Error (RMSE)
  res$rmse <- sqrt(mean((predict(model) - y_train_data)^2))
  # Goodness-of-Fit Statistics (R-squared)
  res$r_squared <- summary(model)$r.squared

  # Make predictions
  # predictions <- stats::predict(model)
  # Model performance
  # res$rmse_test = caret::RMSE(predictions, test.data[,dependent_variable])
  # res$r_squared_test = caret::R2(predictions, test.data[,dependent_variable])

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(model))
  # conf_int <- confint(model)

  # for a and b extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- paste0("EXP_",row_name,"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name

    exp_coef_estimate <- data.frame(exp_a_estimate = coefficients[i,1])
    colnames(exp_coef_estimate) <- paste0("EXP_",row_name,"_ESTIMATE",sep="")

    res <- cbind(res, p_value, exp_coef_estimate)
  }

  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  return (res)
}

# exp_model("exp_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
