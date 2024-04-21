polynomial_model <- function (family_test, tempDataFrame, sig.formula)
{

  # plynomial_degree_partition-partition_percentage
  plynomial_params <- unlist(strsplit(as.character(family_test),"_"))

  degree <- as.numeric(plynomial_params[2])
  res <- data.frame("PL_DEGREE"= degree)
  partition_percentage <- as.numeric(plynomial_params[3])
  res$PL_PERC <- partition_percentage

  if(lengh(plynomial_params)==4)
    res$r_model <- paste0("polynomial",plynomial_params[4], sep="")
  else
    res$r_model <- "polynomial"


  tempDataFrame <- as.data.frame(tempDataFrame)
  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  if(lengh(plynomial_params)==4)
  {
    if(param[4]=="predictor")
    {
      dependent_variable <- dep_var[[3]]
      independent_variable <- dep_var[[2]]
    }
    else
    {
      dependent_variable <- dep_var[[2]]
      independent_variable <- dep_var[[3]]
    }
  }

  # Split the data into training and test set
  training.samples <- tempDataFrame[, dependent_variable] %>% caret::createDataPartition(p = partition_percentage, list = FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # Build the model
  model <- stats::lm(eval(parse(text=dependent_variable)) ~ stats::poly(eval(parse(text=independent_variable)), degree, raw = TRUE), data = train.data)
  # model <- stats::lm(train.data[,dependent_variable] ~ stats::poly(train.data[,independent_variable], degree, raw = TRUE))
  # dependent_variable_train <- train.data[,dependent_variable]
  # independent_variable_train <- train.data[,independent_variable]
  # model <- stats::lm(dependent_variable_train ~ stats::poly(independent_variable_train, degree, raw = TRUE))

  # Calculate Adjusted R-squared
  res$adj_r_squared <- summary(model)$adj.r.squared
  # Calculate Sum of Squares due to Error (SSE)
  res$sse <- sum(model$residuals^2)
  # Calculate Root Mean Square Error (RMSE)
  res$rmse <- sqrt(mean(model$residuals^2))
  # Goodness-of-Fit Statistics (R-squared)
  res$r_squared <- summary(model)$r.squared

  # Make predictions
  predictions <- model %>% stats::predict(test.data)
  # Model performance
  res$rmse_test = caret::RMSE(predictions, test.data[,dependent_variable])
  res$r_squared_test = caret::R2(predictions, test.data[,dependent_variable])

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(model))
  # conf_int <- confint(model)

  # for each group combination extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    if (row_name =="(Intercept)")
      pval_name <- "PL_INTERCEPT_PVALUE"
    else
      pval_name <- paste0("PL_DEGREE_",as.character(i -1 ),"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name
    res <- cbind(res, p_value)
  }

  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  return (res)
}

# polynomial_model("polynomial_4_0.8", sample_sheet, "Age ~ DELTARQ_HYPO")
