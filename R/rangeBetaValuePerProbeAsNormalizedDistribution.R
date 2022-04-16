#--- rangeBetaValuePerProbeAsNormalizedDistribution ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' Nrmalize the methyation level distributing as a normalized distribution
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution

rangeBetaValuePerProbeAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 3) {

  populationMatrixDim <- dim(populationMatrix)
  beta_values <- populationMatrix[, 2:populationMatrixDim[2]]
  row.names(beta_values) <- populationMatrix[, 1]

  betaQ1Values <- future.apply::future_apply( beta_values, 1, stats::quantile, 0.25)
  betaQ3Values <- future.apply::future_apply( beta_values, 1, stats::quantile, 0.75)

  betaMedianValues <- future.apply::future_apply( beta_values, 1, stats::median)
  betaValuesIQR <- future.apply::future_apply( beta_values, 1, stats::IQR)

  if (!test_match_order(row.names(beta_values), row.names(betaMedianValues))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  if (!test_match_order(row.names(beta_values), row.names(betaValuesIQR))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

  row.names(beta_inferior_thresholds) <- row.names(betaMedianValues)

  result <- list(beta_inferior_thresholds = beta_inferior_thresholds, beta_superior_thresholds = beta_superior_thresholds, betaMedianValues = betaMedianValues)

  gc()
  return(result)
}
