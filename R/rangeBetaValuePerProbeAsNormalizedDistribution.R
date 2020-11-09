#--- rangeBetaValuePerProbeAsNormalizedDistribution ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' Nrmalize the methyation level distributing as a normalized distribution
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution
#' @export
#' @examples
#' rangeBetaValuePerProbeAsNormalizedDistribution(
#' populationMatrix = controlPopulationMatrix,
#' iqrTimes = 3
#' )
rangeBetaValuePerProbeAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 3) {


  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1)
  registerDoParallel(computationCluster)

  populationMatrixDim <- dim(populationMatrix)
  betaValues <- populationMatrix[, 2:populationMatrixDim[2]]
  row.names(betaValues) <- populationMatrix[, 1]

  betaQ1Values <- parApply(computationCluster, betaValues, 1, quantile, 0.25)
  betaQ3Values <- parApply(computationCluster, betaValues, 1, quantile, 0.75)

  betaMedianValues <- parApply(computationCluster, betaValues, 1, median)
  betaValuesIQR <- parApply(computationCluster, betaValues, 1, IQR)

  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  if (!test_match_order(row.names(betaValues), row.names(betaValuesIQR)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  betaInferiorThresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)

  result <- list(betaInferiorThresholds = betaInferiorThresholds, betaSuperiorThresholds = (betaQ3Values + (iqrTimes * betaValuesIQR)), betaMedianValues = betaMedianValues)

  stopCluster(computationCluster)
  return(result)
}