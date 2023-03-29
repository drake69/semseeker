#' Title
#'
#' @param family_test
#' @param tempDataFrame
#' @param sig.formula
#'
#' @return
#' @export
#'
#' @examples
glm_model <- function(family_test, tempDataFrame, sig.formula )
{

  result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))

  pvalue <- summary(result_glm )$coeff[-1, 4][1]
  beta_value <- (summary(result_glm )$coeff[-1, 1][1])
  aic_value <- (result_glm$aic)
  residuals <-  result_glm$resid
  #calculate shapiro of working residuals
  shapiro_pvalue <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
  # Breusch_Pagan_pvalue <- lmtest::bptest( data=residuals )$p.value
  ci.lower <- NA
  ci.upper <- NA
  r_model <- "stats::glm"
  std.error <- "NA"
  n_permutations <- NA
  ci.lower.adjusted <- NA
  ci.upper.adjusted <- NA

  return (data.frame(ci.lower,ci.upper, pvalue, beta_value,aic_value,residuals,shapiro_pvalue, r_model,std.error, n_permutations,ci.lower.adjusted,ci.upper.adjusted ))
}
