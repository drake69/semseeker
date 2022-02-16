#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param methylationData whole matrix of data to analyze.
#' @param slidingWindowSize size of the sliding widows to compute epilesions
#' default 11 probes.
#' @param betaSuperiorThresholds  data frame to select, from the sample sheet,
#' samples to use as control as study population and as refereces two vectors
#' within the first vector the names of the selection colum and tha second
#' vector the study population selector,
#' @param betaInferiorThresholds name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param sampleSheet name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param betaMedians name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param bonferroniThreshold threshold to define which pValue accept for
#' @param probeFeatures probes detail from 27 to EPIC illumina dataset
#' lesions definition
#' @return files into the result folder with pivot table and bedgraph.

#' @importFrom foreach %dopar%
analizePopulation <- function(methylationData, slidingWindowSize, betaSuperiorThresholds, betaInferiorThresholds, sampleSheet, betaMedians, bonferroniThreshold = 0.05, probeFeatures) {

  # browser()
  start_time <- Sys.time()
  message("AnalyzePopulation warmingUP ", Sys.time())

  ### get beta_values ########################################################
  sampleSheet <- sampleSheet[order(sampleSheet[, "Sample_ID"], decreasing = FALSE), ]
  existentSamples <- colnames(methylationData)
  sampleNames <- sampleSheet$Sample_ID
  sampleToSelect <- existentSamples[sampleNames %in% existentSamples]
  missedSample <- setdiff(setdiff(sampleNames, existentSamples), "PROBE")

  if (length(missedSample) != 0) {
    message("These samples data are missed: ", paste0(missedSample, sep = " "), Sys.time())
  }

  message("WarmedUP AnalyzePopulation", Sys.time())
  message("Start population analyze ", Sys.time())

  summaryFileName <- file.path(ssEnv$resultFolderData, "summary.csv")

  i <- 0
  # summaryPopulation <-  foreach::foreach( i =1:nrow(sampleSheet), .combine='rbind', .export = ssEnv$functionToExport) %dopar% {
  for(i in 1:nrow(sampleSheet) ) {
    localSampleDetail <- sampleSheet[i,]
    betaValues <- methylationData[, localSampleDetail$Sample_ID]
    hyperResult <- analyzeSingleSample(envir=ssEnv, values = betaValues, slidingWindowSize = slidingWindowSize,  thresholds = betaSuperiorThresholds, figure="HYPER", sampleDetail = localSampleDetail, bonferroniThreshold = bonferroniThreshold, probeFeatures = probeFeatures)
    hypoResult <- analyzeSingleSample(envir=ssEnv, values = betaValues, slidingWindowSize = slidingWindowSize,  thresholds = betaInferiorThresholds, figure="HYPO", sampleDetail = localSampleDetail, bonferroniThreshold = bonferroniThreshold, probeFeatures = probeFeatures)
    bothResult <- analyzeSingleSampleBoth(envir=ssEnv, sampleDetail =  localSampleDetail)
    deltaResult <- deltaSingleSample (envir = ssEnv, values = betaValues, highThresholds = betaSuperiorThresholds, lowThresholds = betaInferiorThresholds, sampleDetail = localSampleDetail, betaMedians = betaMedians, probeFeatures = probeFeatures)
    sampleStatusTemp <- c( "Sample_ID"=localSampleDetail$Sample_ID, deltaResult, hyperResult, hypoResult, bothResult)
    sampleStatusTemp <- data.frame(sampleStatusTemp)
    sampleStatusTemp <- data.frame(t(sampleStatusTemp))
    rownames(sampleStatusTemp) <- c(localSampleDetail$Sample_ID)
    as.data.frame(sampleStatusTemp)
    if(exists("summaryPopulation"))
      summaryPopulation <- rbind(summaryPopulation, sampleStatusTemp)
    else
      summaryPopulation <- sampleStatusTemp
  }

  message("Row count result:", nrow(summaryPopulation))
  rm(methylationData)
  gc()

  message("Completed population analysis ", Sys.time())
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  message("Completed population with Excel summary", Sys.time(), " Time taken: ", time_taken)

  return(summaryPopulation)
}


