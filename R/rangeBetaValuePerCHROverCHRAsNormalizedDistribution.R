

#' Title
#'
#' @param populationMatrix
#' @param iqrTimes
#' @param probeFeatures
#'
#' @return
#' @export
#'
#' @examples
rangeBetaValuePerCHROverCHRAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 1, probeFeatures) {

  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1)
  registerDoParallel(computationCluster)

  betaValues <- populationMatrix[, -(1)]

  betaValuesMergedComplete <- data.frame(gene = probeFeatures$gene, CHR = probeFeatures$CHR, START = probeFeatures$START, PROBE = probeFeatures$Probe, betaValues)
  betaValuesMergedNoCHR <- subset(betaValuesMergedComplete, is.na(betaValuesMergedComplete$CHR) | is.null(betaValuesMergedComplete$CHR) | betaValuesMergedComplete$CHR == "")

  ## by CHR median #
  betaValuesMergedCHR <- subset(betaValuesMergedComplete, (!is.na(betaValuesMergedComplete$CHR) & !is.null(betaValuesMergedComplete$CHR) & betaValuesMergedComplete$CHR != ""))

  # # browser()
  betaValuesNoCHRMatrix <- betaValuesMergedNoCHR[, -(1:4)]

  betaMedianValuesNoCHR <- parApply(computationCluster, betaValuesNoCHRMatrix, 1, median)
  betaMedianValuesNoCHR <- data.frame(betaMedianValuesNoCHR)
  colnames(betaMedianValuesNoCHR) <- "median"

  betaIQRValuesNoCHR <- parApply(computationCluster, betaValuesNoCHRMatrix, 1, IQR)
  betaIQRValuesNoCHR <- data.frame(betaIQRValuesNoCHR)
  colnames(betaIQRValuesNoCHR) <- "iqr"

  ## browser()

  chromosomes <- as.character(unique(attributes(betaValuesMergedCHR$CHR)$levels))
  chromosomes <- subset(chromosomes, chromosomes != "")
  # betaMedianIQRValues <- foreach(i = 1:length(chromosomes), .combine = rbind) %dopar%
  for (i in 1:length(chromosomes)) {
    geneSelected <- chromosomes[i]
    betaValuesCHR <- subset(betaValuesMergedCHR, CHR == geneSelected)
    probes <- as.character(betaValuesCHR[, "PROBE"])
    betaValuesCHRAsMatrix <- as.matrix(betaValuesCHR[, -(1:4)], ncol = 1)
    betaValuesCHRAsMatrix <- unmatrix(betaValuesCHRAsMatrix)
    # print(e1071::skewness(betaValuesCHRAsMatrix)) hist(betaValuesCHRAsMatrix) fw <- fitdistrplus::descdist(betaValuesCHRAsMatrix) print(summary(fw))
    medianValues <- rep(median(betaValuesCHRAsMatrix), dim(betaValuesCHR)[1])
    iqrValues <- rep(IQR(betaValuesCHRAsMatrix), dim(betaValuesCHR)[1])
    data.frame(row.names = probes, median = medianValues, iqr = iqrValues)
  }

  ## browser()

  betaMedianValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[, "median"])
  colnames(betaMedianValues) <- "median"

  betaMedianValues <- rbind(betaMedianValues, betaMedianValuesNoCHR)
  betaMedianValues <- betaMedianValues[order(row.names(betaMedianValues), decreasing = FALSE), ]
  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  betaIQRValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[, "iqr"])
  colnames(betaIQRValues) <- "iqr"

  betaIQRValues <- rbind(betaIQRValues, betaIQRValuesNoCHR)
  betaIQRValues <- betaIQRValues[order(row.names(betaIQRValues), decreasing = FALSE), ]
  if (!test_match_order(row.names(betaValues), row.names(betaIQRValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  # # browser()

  betaInferiorThresholds <- (betaMedianValues - (iqrTimes * betaIQRValues))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)

  result <- list(betaInferiorThresholds = betaInferiorThresholds, betaSuperiorThresholds = (betaMedianValues + (iqrTimes * betaIQRValues)), betaMedianValues = betaMedianValues)

  stopCluster(computationCluster)
  return(result)
}
