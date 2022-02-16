#--- rangeBetaValuePerProbeAsNormalizedDistribution ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' Nrmalize the methyation level distributing as a normalized distribution
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution

rangeBetaValuePerProbeAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 3) {

  populationMatrixDim <- dim(populationMatrix)
  betaValues <- populationMatrix[, 2:populationMatrixDim[2]]
  row.names(betaValues) <- populationMatrix[, 1]

  betaQ1Values <- future.apply::future_apply( betaValues, 1, stats::quantile, 0.25)
  betaQ3Values <- future.apply::future_apply( betaValues, 1, stats::quantile, 0.75)

  betaMedianValues <- future.apply::future_apply( betaValues, 1, stats::median)
  betaValuesIQR <- future.apply::future_apply( betaValues, 1, stats::IQR)

  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  if (!test_match_order(row.names(betaValues), row.names(betaValuesIQR))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  betaInferiorThresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  betaSuperiorThresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)

  result <- list(betaInferiorThresholds = betaInferiorThresholds, betaSuperiorThresholds = betaSuperiorThresholds, betaMedianValues = betaMedianValues)

  gc()
  return(result)
}
