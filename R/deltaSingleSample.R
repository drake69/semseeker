#' Title
#'
#' @param values values of methylation
#' @param resultFolder folder to store result files
#' @param highThresholds highest threshold to use for comparison
#' @param lowThresholds lowest threshold to use for comparison
#' @param sampleName name of sample to analyze
#' @param betaMedians median to use for calculation
#' @param subFileExtension extension to pre pend to file name
#'
#' @return none

#'
deltaSingleSample <- function(values, resultFolder, highThresholds, lowThresholds, sampleName, betaMedians, subFileExtension, probeFeatures) {

  MUTATIONS <- NULL
  message(sampleName, " ", "... Sample analysis warmingUP ", Sys.time())

  colnames(values) <- "VALUE"

  ### get probeFeatures ################################################################################################

  message(sampleName, " ", "Sample analysis WarmedUP ...", Sys.time())
  message(sampleName, " ", "Start sample analyze ", Sys.time())

  ### get probesOverThreshold ################################################################################################

  mutationAbove <- values > highThresholds
  mutationBelow <- values < lowThresholds
  mutation <- (mutationBelow + mutationAbove) > 0
  colnames(mutation) <- "MUTATIONS"

  message(sampleName, " ", "Got outliers ", Sys.time())

  ### get deltas #########################################################

  deltas <- values - betaMedians
  colnames(deltas) <- "DELTA"

  if (!test_match_order(row.names(mutation), probeFeatures$PROBE)) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }
  if (!test_match_order(row.names(deltas), probeFeatures$PROBE)) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  deltasAnnotated <- data.frame(as.data.frame(probeFeatures), deltas, "MUTATIONS" = mutation)

  deltasAnnotatedSorted <- sortByCHRandSTART(deltasAnnotated)

  deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, MUTATIONS == 1)[, c("CHR", "START", "END", "DELTA")]

  dumpSampleAsBedFile(
    dataToDump = deltasAnnotatedSorted,
    fileExtension = ".DELTAS.METHYLATION.bedgraph",
    resultFolder = resultFolder,
    resultSubFolder = "DELTAS_METHYLATION",
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "DELTA", "SAMPLENAME")
  )
}
