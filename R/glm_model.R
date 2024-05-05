#' Title
#'
#' @param family_test regression model to apply
#' @param tempDataFrame data frame to use for the model
#' @param sig.formula formula to apply the model
#'
glm_model <- function(family_test, tempDataFrame, sig.formula)
{
  result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
  res <- data.frame("pvalue" = summary(result_glm)$coeff[-1, 4][1])
  # browser()
  res$statistic_parameter <- (summary(result_glm)$coeff[-1, 1][1])
  res$aic_value <- (result_glm$aic)
  res$std.error <- summary(result_glm)$coeff[-1, 2][1]

  # box.plot(dataFrameToPlot, independent_variable,burdenValue, transformation, family_test)
  # calculate rmse

  if(family_test=="gaussian")
  {
    residuals <-  result_glm$resid
    res$rmse <- sqrt(mean(result_glm$residuals^2))
    #calculate shapiro of working residuals
    res$shapiro_pvalue_residuals <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
    # calculate residuals homoschedadasticity
    res$residuals_breusch_pagan_pvalue_compl <- 1- lmtest::bptest(result_glm)$p.value
    # bartlett of residuals
    # res$residuals_bartlett_pvalue_compl <- ifelse((length(residuals)>3 & length(unique(residuals))>3),1-(stats::bartlett.test(residuals)$p.value),NA)
  }

  #  calculate r-squared
  # res$r_squared <- summary(result_glm)$r.squared

  res$r_model <- "stats::glm"
  return (res)
}
