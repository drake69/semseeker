#' Title
#'
#' @param family_test regression model to apply
#' @param tempDataFrame data frame to use for the model
#' @param sig.formula formula to apply the model
#'
glm_model <- function(family_test, tempDataFrame, sig.formula )
{
  result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
  res <- data.frame("pvalue" = summary(result_glm)$coeff[-1, 4][1])
  res$statistic_parameter <- (summary(result_glm)$coeff[-1, 1][1])
  res$aic_value <- (result_glm$aic)
  if(family_test=="gaussian")
  {
    residuals <-  result_glm$resid
    #calculate shapiro of working residuals
    res$shapiro_pvalue_residuals <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
    res$Breusch_Pagan_pvalue_residuals <- lmtest::bptest( data=residuals )$p.value
    # bartlett of residuals
    res$bartlett_pvalue_residuals <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::bartlett.test(residuals, as.factor(tempDataFrame[, independent_variable]))$p.value) else NA
  }
  res$r_model <- "stats::glm"
  return (res)
}
