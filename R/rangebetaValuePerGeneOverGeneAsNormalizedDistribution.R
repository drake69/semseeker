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
rangeBetaValuePerGeneOverGeneAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 1, probeFeatures) {

  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1)
  registerDoParallel(computationCluster)

  betaValues <- populationMatrix[, -(1)]

  betaValuesMergedComplete <- data.frame(CHR = probeFeatures$CHR, gene = probeFeatures$gene, START = probeFeatures$START, PROBE = probeFeatures$Probe, betaValues)
  betaValuesMergedNoGene <- subset(betaValuesMergedComplete, is.na(betaValuesMergedComplete$gene) | is.null(betaValuesMergedComplete$gene) | betaValuesMergedComplete$gene == "")

  ## by gene median #
  betaValuesMergedGene <- subset(betaValuesMergedComplete, (!is.na(betaValuesMergedComplete$gene) & !is.null(betaValuesMergedComplete$gene) & betaValuesMergedComplete$gene != ""))


  # # browser()
  betaValuesNoGeneMatrix <- betaValuesMergedNoGene[, -(1:4)]

  betaMedianValuesNoGene <- parApply(computationCluster, betaValuesNoGeneMatrix, 1, median)
  betaMedianValuesNoGene <- data.frame(betaMedianValuesNoGene)
  colnames(betaMedianValuesNoGene) <- "median"

  betaIQRValuesNoGene <- parApply(computationCluster, betaValuesNoGeneMatrix, 1, IQR)
  betaIQRValuesNoGene <- data.frame(betaIQRValuesNoGene)
  colnames(betaIQRValuesNoGene) <- "iqr"

  # # browser()

  genes <- as.character(unique(attributes(betaValuesMergedGene$gene)$levels))
  genes <- subset(genes, genes != "")
  betaMedianIQRValues <- foreach(i = 1:length(genes), .combine = rbind) %dopar% # for (i in 1:length(genes))
    {
      geneSelected <- genes[i]
      betaValuesGene <- subset(betaValuesMergedGene, gene == geneSelected)
      probes <- as.character(betaValuesGene[, "PROBE"])
      betaValuesGeneAsMatrix <- as.matrix(betaValuesGene[, -(1:4)])

      medianValues <- rep(median(betaValuesGeneAsMatrix), dim(betaValuesGene)[1])
      iqrValues <- rep(IQR(betaValuesGeneAsMatrix), dim(betaValuesGene)[1])
      data.frame(row.names = probes, median = medianValues, iqr = iqrValues)
    }

  # # browser()

  betaMedianValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[, "median"])
  colnames(betaMedianValues) <- "median"

  betaMedianValues <- rbind(betaMedianValues, betaMedianValuesNoGene)
  betaMedianValues <- betaMedianValues[order(row.names(betaMedianValues), decreasing = FALSE), ]
  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  betaIQRValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[, "iqr"])
  colnames(betaIQRValues) <- "iqr"

  betaIQRValues <- rbind(betaIQRValues, betaIQRValuesNoGene)
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
