#' Title
#'
#' @param family_test
#' @param tempDataFrame
#' @param sig.formula
#' @param burdenValue
#' @param independent_variable
#'
#' @return
#' @export
#'
#' @examples
test_model <- function (family_test, tempDataFrame, sig.formula,burdenValue,independent_variable )
{
  if(family_test=="wilcoxon")
  {
    result_w  <- suppressWarnings(stats::wilcox.test(formula= sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE))
    pvalue <- result_w$p.value
    r_model <- "stats_wilcox.test"
  }

  if(family_test=="t.test")
  {
    result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame))
    pvalue <- result_w$p.value
    r_model <- "stats_t.test"
  }

  if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
  {
    result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
    pvalue <- result_cor$p.value
    r_model <- "stats_cor.test"
    beta_value <- result_cor$estimate
  }

  ci.lower <- NA
  ci.upper <- NA
  beta_value <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  n_permutations <- NA
  ci.lower.adjusted <- NA
  ci.upper.adjusted <- NA

  return (data.frame(ci.lower,ci.upper, pvalue, beta_value,aic_value,residuals,shapiro_pvalue, r_model,std.error,n_permutations,ci.lower.adjusted,ci.upper.adjusted))

}
