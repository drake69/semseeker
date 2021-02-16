#--- rangeBetaValuePerProbeAsNormalizedDistribution ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' Nrmalize the methyation level distributing as a normalized distribution
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution

rangeBetaValuePerProbeAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 3) {
  computationCluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1)
  doParallel::registerDoParallel(computationCluster)

  populationMatrixDim <- dim(populationMatrix)
  betaValues <- populationMatrix[, 2:populationMatrixDim[2]]
  row.names(betaValues) <- populationMatrix[, 1]

  betaQ1Values <- parallel::parApply(computationCluster, betaValues, 1, stats::quantile, 0.25)
  betaQ3Values <- parallel::parApply(computationCluster, betaValues, 1, stats::quantile, 0.75)

  betaMedianValues <- parallel::parApply(computationCluster, betaValues, 1, stats::median)
  betaValuesIQR <- parallel::parApply(computationCluster, betaValues, 1, stats::IQR)

  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  if (!test_match_order(row.names(betaValues), row.names(betaValuesIQR))) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  betaInferiorThresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)

  result <- list(betaInferiorThresholds = betaInferiorThresholds, betaSuperiorThresholds = (betaQ3Values + (iqrTimes * betaValuesIQR)), betaMedianValues = betaMedianValues)

  parallel::stopCluster(computationCluster)
  return(result)
}
