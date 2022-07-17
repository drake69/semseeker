#--- range_beta_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of beta values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution

range_beta_values <- function(populationMatrix, iqrTimes = 3) {

  populationMatrixDim <- dim(populationMatrix)
  beta_values <- as.data.frame(populationMatrix[, 2:populationMatrixDim[2]])
  row.names(beta_values) <- populationMatrix[, 1]

  betaQ1Values <- as.data.frame(future.apply::future_apply( beta_values, 1, stats::quantile, 0.25))
  betaQ3Values <- as.data.frame(future.apply::future_apply( beta_values, 1, stats::quantile, 0.75))

  beta_median_values <- as.data.frame(future.apply::future_apply( beta_values, 1, stats::median))
  colnames(beta_median_values) <- c("beta_median_values")

  betaValuesIQR <- as.data.frame(future.apply::future_apply( beta_values, 1, stats::IQR))

  if (!test_match_order(row.names(beta_values), row.names(beta_median_values))) {
    stop("Wrong order matching Probes and BetaMedianvalues!", Sys.time())
  }

  if (!test_match_order(row.names(beta_values), row.names(betaValuesIQR))) {
    stop("Wrong order matching Probes and BetaMedianvalues IQR!", Sys.time())
  }

  beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

  row.names(beta_inferior_thresholds) <- row.names(beta_median_values)
  colnames(beta_inferior_thresholds) <- c("beta_inferior_thresholds")

  row.names(beta_superior_thresholds) <- row.names(beta_median_values)
  colnames(beta_superior_thresholds) <- c("beta_superior_thresholds")

  result <- list(beta_inferior_thresholds = beta_inferior_thresholds, beta_superior_thresholds = beta_superior_thresholds, beta_median_values = beta_median_values)

  gc()
  return(result)
}
