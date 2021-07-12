#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param methylationData whole matrix of data to analyze.
#' @param slidingWindowSize size of the sliding widows to compute epilesions
#' default 11 probes.
#' @param resultFolder folder to save computed epimutations and bedgraphs files.
#' @param logFolder folder to output log
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
#' @param populationName name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param bonferroniThreshold threshold to define which pValue accept for
#' @param probeFeatures probes detail from 27 to EPIC illumina dataset
#' lesions definition
#' @return files into the result folder with pivot table and bedgraph.

#' @importFrom foreach %dopar%
analizePopulation <- function(methylationData, slidingWindowSize, resultFolder, logFolder, betaSuperiorThresholds, betaInferiorThresholds, sampleSheet, betaMedians, populationName, bonferroniThreshold = 0.05, probeFeatures) {

  if (.Platform$OS.type == "windows") {
    withAutoprint({
      utils::memory.size()
      utils::memory.size(TRUE)
      utils::memory.limit(16000)
    })
  }

  start_time <- Sys.time()
  message("... warmingUP ", Sys.time())

  if (logFolder != "" && !dir.exists(logFolder)) {
    dir.create(logFolder)
  }

  resultFolderBase <- resultFolder
  # population folder check
  resultFolder <- paste0(resultFolder, "/", populationName)
  if (resultFolder != "" && !dir.exists(resultFolder)) {
    dir.create(resultFolder)
  }

  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
  outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
  print(outFile)
  computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computation_cluster)

  # options(digits = 22)
  parallel::clusterExport(envir=environment(), cl = computation_cluster,
                          varlist = list( "analyzeSingleSample", "dumpSampleAsBedFile", "deltaSingleSample",
                                      "createPivotResultFromMultipleBed", "sortByCHRandSTART", "test_match_order", "getLesions", "addCellToDataFrame",
                                      "PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body"))

  ### get beta_values ########################################################
  sampleSheet <- sampleSheet[order(sampleSheet[, "Sample_ID"], decreasing = FALSE), ]
  sampleSheet[, "Probes_Count"] <- 0
  sampleSheet[, "Hyper_Mutations"] <- 0
  sampleSheet[, "Hyper_Lesions"] <- 0
  sampleSheet[, "Hypo_Mutations"] <- 0
  sampleSheet[, "Hypo_Lesions"] <- 0

  existentSamples <- colnames(methylationData)
  sampleNames <- sampleSheet$Sample_ID
  sampleToSelect <- existentSamples[sampleNames %in% existentSamples]
  missedSample <- setdiff(setdiff(sampleNames, existentSamples), "PROBE")

  if (length(missedSample) != 0) {
    message("These samples data are missed: ", paste0(missedSample, sep = " "), Sys.time())
  }

  message("WarmedUP ...", Sys.time())
  message("Start population analyze ", Sys.time())

  summaryFileName <- file.path(resultFolder, "summary.csv")
  system(paste0("echo '", paste(colnames(sampleSheet), collapse = ","), "' > ", summaryFileName, sep = ""))

  columnIndexes <- which((colnames(methylationData)%in%sampleSheet$Sample_ID))

  # browser()
  i <- NULL
  # for (i in columnIndexes) {


  splittedColumnIndexes <- split(columnIndexes, ceiling(seq_along(columnIndexes)/200))

  # availableRam <- memuse::Sys.meminfo()$freeram@size
  # mDataSize <-  pryr::object_size(methylationData)
  # maxDataSize <- availableRam * 0.9 /nCore

  for (x in 1:length(splittedColumnIndexes))
  {
    # browser()
    betaValues <- methylationData[,as.vector(unlist(splittedColumnIndexes[x]))]

    # betaValuesHyper <- betaValues[betaValues > betaSuperiorThresholds]  - betaSuperiorThresholds
    # betaValuesHypo <- betaInferiorThresholds - betaValues[betaValues < betaInferiorThresholds]

    foreach::foreach(i = 1:dim(betaValues)[2] , .packages=c("dplyr")) %dopar% {
    # for(i in 1:dim(betaValues)[2] ) {

      sampleName <- colnames(betaValues)[i]
      # browser()
      message("Starting sample analysis number: ", i, " ", Sys.time())
      sampleDetail <- sampleSheet[sampleSheet$Sample_ID==sampleName, ]
      colnames(sampleDetail) <- colnames(sampleSheet)

      deltaSingleSample(
        probeFeatures = probeFeatures,  values = betaValues[i], resultFolder = resultFolder, highThresholds = betaSuperiorThresholds, lowThresholds = betaInferiorThresholds, sampleName = sampleDetail$Sample_ID,
        betaMedians = betaMedians, subFileExtension = "DELTAS"
      )

      sampleStatusTemp <- analyzeSingleSample(
        values = betaValues[i], slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaSuperiorThresholds, comparison = `>`, sampleName = sampleDetail$Sample_ID,
        subFileExtension = "HYPER", bonferroniThreshold = bonferroniThreshold, probeFeatures = probeFeatures
      )

      # browser()
      sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Probes_Count", cellValue = sampleStatusTemp["probesCount"])
      sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hyper_Mutations", cellValue = sampleStatusTemp["mutationCount"])
      sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hyper_Lesions", cellValue = sampleStatusTemp["lesionCount"])

      sampleStatusTemp <- analyzeSingleSample(
        values = betaValues[i], slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaInferiorThresholds, comparison = `<`, sampleName = sampleDetail$Sample_ID,
        subFileExtension = "HYPO", probeFeatures = probeFeatures
      )

      # # browser()
      sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hypo_Mutations", cellValue = sampleStatusTemp["mutationCount"])
      sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hypo_Lesions", cellValue = sampleStatusTemp["lesionCount"])

      tempExtension <- stringi::stri_rand_strings(1, 10)
      filePath <- paste0(summaryFileName, tempExtension, sep = "")
      utils::write.table(sampleDetail, file = filePath, quote = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

      if (.Platform$OS.type == "windows") {

        command <- paste0("type ", (filePath), " > ",(summaryFileName), sep = "")
        command <- gsub ("/","\\\\",command)
        shell(command, intern = TRUE)

        command <- paste0("del ", (filePath), sep = "")
        command <- gsub ("/","\\\\",command)
        #print(command)
        shell(command, intern = TRUE)
        # system2(paste0("type ", shQuote(filePath), " > ",shQuote(summaryFileName), sep = ""))
        # system2(paste0("rm ", filePath, sep = ""))
      } else
      {
        system(paste0("cat ", filePath, " >> ", summaryFileName, sep = ""))
        system(paste0("rm ", filePath, sep = ""))
      }

      gc()
    }

    rm(betaValues)
    gc()
  }


  # myTest <- function(betaData)
  # {
  #     sampleName <- colnames(betadata)
  #     # browser()
  #     message("Starting sample analysis number: ", i, " ", Sys.time())
  #     sampleDetail <- sampleSheet[sampleSheet$Sample_ID==sampleName, ]
  #     colnames(sampleDetail) <- colnames(sampleSheet)
  #
  #     deltaSingleSample(
  #       probeFeatures = probeFeatures,  values = betaData, resultFolder = resultFolder, highThresholds = betaSuperiorThresholds, lowThresholds = betaInferiorThresholds, sampleName = sampleDetail$Sample_ID,
  #       betaMedians = betaMedians, subFileExtension = "DELTAS"
  #     )
  #
  #     sampleStatusTemp <- analyzeSingleSample(
  #       values = betaData, slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaSuperiorThresholds, comparison = `>`, sampleName = sampleDetail$Sample_ID,
  #       subFileExtension = "HYPER", bonferroniThreshold = bonferroniThreshold, probeFeatures = probeFeatures
  #     )
  #
  #     # browser()
  #     sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Probes_Count", cellValue = sampleStatusTemp["probesCount"])
  #     sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hyper_Mutations", cellValue = sampleStatusTemp["mutationCount"])
  #     sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hyper_Lesions", cellValue = sampleStatusTemp["lesionCount"])
  #
  #     sampleStatusTemp <- analyzeSingleSample(
  #       values = betaData, slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaInferiorThresholds, comparison = `<`, sampleName = sampleDetail$Sample_ID,
  #       subFileExtension = "HYPO", probeFeatures = probeFeatures
  #     )
  #
  #     # # browser()
  #     sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hypo_Mutations", cellValue = sampleStatusTemp["mutationCount"])
  #     sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_ID", cellValueSelection = sampleDetail$Sample_ID, colname = "Hypo_Lesions", cellValue = sampleStatusTemp["lesionCount"])
  #
  #     tempExtension <- stringi::stri_rand_strings(1, 10)
  #     filePath <- paste0(summaryFileName, tempExtension, sep = "")
  #     utils::write.table(sampleDetail, file = filePath, quote = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  #
  #     if (.Platform$OS.type == "windows") {
  #
  #       command <- paste0("type ", (filePath), " > ",(summaryFileName), sep = "")
  #       command <- gsub ("/","\\\\",command)
  #       shell(command, intern = TRUE)
  #
  #       command <- paste0("del ", (filePath), sep = "")
  #       command <- gsub ("/","\\\\",command)
  #       #print(command)
  #       shell(command, intern = TRUE)
  #       # system2(paste0("type ", shQuote(filePath), " > ",shQuote(summaryFileName), sep = ""))
  #       # system2(paste0("rm ", filePath, sep = ""))
  #     } else
  #     {
  #       system(paste0("cat ", filePath, " >> ", summaryFileName, sep = ""))
  #       system(paste0("rm ", filePath, sep = ""))
  #     }
  # }
  #
  #
  # parallel::parApply(computation_cluster, methylationData[, columnIndexes] , 2, myTest)


  sampleSheet <- utils::read.csv(file = summaryFileName)
  #file.remove(summaryFileName)

  ### get countProbesEpiMutatedPerSample ######################################################### probes_above_high_thresholds <- beta_values > betaSuperiorThresholds probes_below_low_thresholds <- beta_values <
  ### betaInferiorThresholds count_probes_above_high_threshold_per_sample <- colSums(probes_above_high_thresholds) count_probes_below_low_threshold_per_sample <- colSums(probes_below_low_thresholds)
  ### count_probes_out_of_range_per_sample <- count_probes_above_high_threshold_per_sample + count_probes_below_low_threshold_per_sample count_probes_epi_mutated_per_sample <- data.frame( 'Sample' = sampleSheet$Sample_ID,
  ### count_probes_below_low_threshold_per_sample, count_probes_above_high_threshold_per_sample, count_probes_out_of_range_per_sample) write.table(count_probes_epi_mutated_per_sample, paste(resultFolder, '/', 'TAB_EPIMUT.txt', sep =
  ### ''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE) message('Got countProbesEpiMutatedPerSample ', Sys.time() ) rm(count_probes_above_high_threshold_per_sample) rm(count_probes_below_low_threshold_per_sample)
  ### rm(count_probes_out_of_range_per_sample) rm(count_probes_epi_mutated_per_sample)

  ### get epiMutationBurden ##############################################################################

  # probesOutOfRange <- (probes_above_high_thresholds) + (probes_below_low_thresholds) epiMutationPerProbes <- rowSums(probesOutOfRange) epiMutationBurden <- data.frame( row.names = probePositions$PROBE, epiMutationPerProbes,
  # probePositions$PROBE) colnames(epiMutationBurden) <- c('Burden', 'PROBE') epiMutationAbovePerProbes <- rowSums(probes_above_high_thresholds) epiMutationAboveBurden <- data.frame(row.names = probePositions$PROBE,
  # epiMutationAbovePerProbes, probePositions$PROBE) colnames(epiMutationAboveBurden) <- c('Burden', 'PROBE') epiMutationBelowBurden <- data.frame(rowSums(probes_below_low_thresholds), probePositions$PROBE)
  # colnames(epiMutationBelowBurden) <- c('Burden', 'PROBE') epiMutation <- data.frame(epiMutationBurden, probesOutOfRange) rm(epiMutationBurden) rm(probesOutOfRange) rm(epiMutationPerProbes) rm(epiMutationAbovePerProbes)
  # message('Got epiMutationBelowBurden ', Sys.time() ) ### get epiMutation ############################################################################## # epiMutationAbove <- data.frame(epiMutationAboveBurden,
  # probes_above_high_thresholds) # epiMutationBelow <- data.frame(epiMutationBelowBurden, probes_below_low_thresholds) epiMutation <- subset(epiMutation, epiMutation$Burden > 0) write.table(epiMutation, paste(resultFolder, '/',
  # 'EPIMUTATIONs.txt', sep = ''), sep = '\t', row.names = FALSE, quote = FALSE) rm(epiMutation) rm(epiMutationBelowBurden) message('Got epiMutation ', Sys.time() )

  parallel::stopCluster(computation_cluster)

  message("Completed population analysis ", Sys.time())

  # xPlot <- colnames(dataToGraph)[-c('SAMPLENAME')] plot( xPlot, dataToGraph[2,2:dim(dataToGraph)[2]] ,type = 'l',col = 'red') for (z in 3:dim(dataToGraph)[1]) { lines(xPlot, dataToGraph[z,2:dim(dataToGraph)[2]], col = 'green') }
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  message("Completed population with Excel summary", Sys.time(), " Time taken: ", time_taken)
}

