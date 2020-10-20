#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param populationMatrix size of the sliding widows to compute epilesions
#' default 11 probes.
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
#' @param probeFeatures name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param betaMeans name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param populationName name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param bonferroniThreshold threshold to define which pValue accept for
#' lesions definition
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#' @examples
#' calculate(
#' populationMatrix
#' slidingWindowSize = 11,
#' resultFolder = resultFolderValue,
#' logFolder = logFolderValue,
#' betaSuperiorThresholds  = ,
#' betaInferiorThresholds =,
#' sampleSheet,
#' probeFeatures,
#' betaMeans,
#' populationName,
#' bonferroniThreshold = 0.05
#' )
analizePopulation <- function(populationMatrix, slidingWindowSize, resultFolder, logFolder, betaSuperiorThresholds, betaInferiorThresholds, sampleSheet, probeFeatures, betaMeans, populationName, bonferroniThreshold = 0.05) {

    if (.Platform$OS.type == "windows")
        withAutoprint({
            memory.size()
            memory.size(TRUE)
            memory.limit(16000)
        })

    start_time <- Sys.time()
    message("... warmingUP ", Sys.time())

    if (logFolder != "" && !dir.exists(logFolder)) {
        dir.create(logFolder)
    }

    # population folder check
    resultFolder <- paste0(resultFolder, "/", populationName)
    if (resultFolder != "" && !dir.exists(resultFolder)) {
        dir.create(resultFolder)
    }

    outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
    computation_cluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
    parallel::registerDoParallel(computation_cluster)

    options(digits = 22)
    clusterExport(computation_cluster, list("analyzeSingleSample", "dumpSampleAsBedFile", "deltaSingleSample", "createPivotResultFromMultipleBed", "sortByCHRandSTART", "mixedsort", "test_match_order", "getLesions", "addCellToDataFrame"))
    ### get probeFeatures ######################################################

    probeFeaturesOriginal <- probeFeatures

    # browser()
    probePositions <- data.frame(PROBE = probeFeatures$Probe, CHR = probeFeatures$CHR, START = probeFeatures$START)
    dim(probePositions)

    ### get beta_values ########################################################
    sampleSheet <- sampleSheet[order(sampleSheet[, "Sample_Name"], decreasing = FALSE), ]
    sampleSheet[, "Probes_Count"] <- 0
    sampleSheet[, "Hyper_Mutations"] <- 0
    sampleSheet[, "Hyper_Lesions"] <- 0
    sampleSheet[, "Hypo_Mutations"] <- 0
    sampleSheet[, "Hypo_Lesions"] <- 0

    existentSamples <- colnames(populationMatrix)
    sampleNames <- sampleSheet$Sample_Name
    sampleToSelect <- existentSamples[existentSamples %in% sampleNames]
    missedSample <- setdiff(setdiff(sampleNames, existentSamples), "PROBE")
    beta_values <- populationMatrix[, sampleToSelect]

    if (length(missedSample) != 0)
        message("These samples data are missed: ", paste0(missedSample, sep = " "), Sys.time())

    message("WarmedUP ...", Sys.time())
    message("Start population analyze ", Sys.time())

    # browser()
    summaryFileName <- paste0(resultFolder, "/", "summary.csv", sep = "")
    system(paste0("echo '", paste(colnames(sampleSheet), collapse = ","), "' > ", summaryFileName, sep = ""))

    # for (i in 1:dim(beta_values)[2]) {
    foreach(i = 1:dim(beta_values)[2]) %dopar% {

        # browser()
        message("Starting sample analysis number: ", i, " ", Sys.time())
        sampleDetail <- sampleSheet[i, ]
        colnames(sampleDetail) <- colnames(sampleSheet)

        deltaSingleSample(values = beta_values[i], resultFolder = resultFolder, highThresholds = betaSuperiorThresholds, lowThresholds = betaInferiorThresholds, sampleName = sampleDetail$Sample_Name, probeFeatures = probeFeatures,
            means = betaMeans, subFileExtension = "DELTAS")

        sampleStatusTemp <- analyzeSingleSample(values = beta_values[i], slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaSuperiorThresholds, comparison = `>`, sampleName = sampleDetail$Sample_Name,
            probeFeatures = probeFeatures, means = betaMeans, subFileExtension = "HYPER", bonferroniThreshold = bonferroniThreshold)

        # browser()
        sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_Name", cellValueSelection = sampleDetail$Sample_Name, colname = "Probes_Count", cellValue = sampleStatusTemp["probesCount"])
        sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_Name", cellValueSelection = sampleDetail$Sample_Name, colname = "Hyper_Mutations", cellValue = sampleStatusTemp["mutationCount"])
        sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_Name", cellValueSelection = sampleDetail$Sample_Name, colname = "Hyper_Lesions", cellValue = sampleStatusTemp["lesionCount"])

        sampleStatusTemp <- analyzeSingleSample(values = beta_values[i], slidingWindowSize = slidingWindowSize, resultFolder = resultFolder, thresholds = betaInferiorThresholds, comparison = `<`, sampleName = sampleDetail$Sample_Name,
            probeFeatures = probeFeatures, means = betaMeans, subFileExtension = "HYPO")

        # browser()
        sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_Name", cellValueSelection = sampleDetail$Sample_Name, colname = "Hypo_Mutations", cellValue = sampleStatusTemp["mutationCount"])
        sampleDetail <- addCellToDataFrame(sampleDetail, colSelection = "Sample_Name", cellValueSelection = sampleDetail$Sample_Name, colname = "Hypo_Lesions", cellValue = sampleStatusTemp["lesionCount"])

        tempExtension <- stringi::stri_rand_strings(1, 10)
        filePath <- paste0(summaryFileName, tempExtension, sep = "")
        write.table(sampleDetail, file = filePath, quote = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
        system(paste0("cat ", filePath, " >> ", summaryFileName, sep = ""))
        system(paste0("rm ", filePath, sep = ""))

        # forEachSampleSheet[i,] <- sampleDetail system(paste0('echo '', paste(sampleDetail,collapse = ',') , '' >> ', summaryFileName, sep = ''))
    }

    sampleSheet <- read.csv(file = summaryFileName)
    # file.remove(summaryFileName)

    ### get countProbesEpiMutatedPerSample ######################################################### probes_above_high_thresholds <- beta_values > betaSuperiorThresholds probes_below_low_thresholds <- beta_values <
    ### betaInferiorThresholds count_probes_above_high_threshold_per_sample <- colSums(probes_above_high_thresholds) count_probes_below_low_threshold_per_sample <- colSums(probes_below_low_thresholds)
    ### count_probes_out_of_range_per_sample <- count_probes_above_high_threshold_per_sample + count_probes_below_low_threshold_per_sample count_probes_epi_mutated_per_sample <- data.frame( 'Sample' = sampleSheet$Sample_Name,
    ### count_probes_below_low_threshold_per_sample, count_probes_above_high_threshold_per_sample, count_probes_out_of_range_per_sample) write.table(count_probes_epi_mutated_per_sample, paste(resultFolder, '/', 'TAB_EPIMUT.txt', sep =
    ### ''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE) message('Got countProbesEpiMutatedPerSample ', Sys.time() ) rm(count_probes_above_high_threshold_per_sample) rm(count_probes_below_low_threshold_per_sample)
    ### rm(count_probes_out_of_range_per_sample) rm(count_probes_epi_mutated_per_sample)

    ### get epiMutationBurden ##############################################################################

    # probesOutOfRange <- (probes_above_high_thresholds) + (probes_below_low_thresholds) epiMutationPerProbes <- rowSums(probesOutOfRange) epiMutationBurden <- data.frame( row.names = probePositions$Probe, epiMutationPerProbes,
    # probePositions$Probe) colnames(epiMutationBurden) <- c('Burden', 'Probe') epiMutationAbovePerProbes <- rowSums(probes_above_high_thresholds) epiMutationAboveBurden <- data.frame(row.names = probePositions$Probe,
    # epiMutationAbovePerProbes, probePositions$Probe) colnames(epiMutationAboveBurden) <- c('Burden', 'Probe') epiMutationBelowBurden <- data.frame(rowSums(probes_below_low_thresholds), probePositions$Probe)
    # colnames(epiMutationBelowBurden) <- c('Burden', 'Probe') epiMutation <- data.frame(epiMutationBurden, probesOutOfRange) rm(epiMutationBurden) rm(probesOutOfRange) rm(epiMutationPerProbes) rm(epiMutationAbovePerProbes)
    # message('Got epiMutationBelowBurden ', Sys.time() ) ### get epiMutation ############################################################################## # epiMutationAbove <- data.frame(epiMutationAboveBurden,
    # probes_above_high_thresholds) # epiMutationBelow <- data.frame(epiMutationBelowBurden, probes_below_low_thresholds) epiMutation <- subset(epiMutation, epiMutation$Burden > 0) write.table(epiMutation, paste(resultFolder, '/',
    # 'EPIMUTATIONs.txt', sep = ''), sep = '\t', row.names = FALSE, quote = FALSE) rm(epiMutation) rm(epiMutationBelowBurden) message('Got epiMutation ', Sys.time() )

    stopCluster(computation_cluster)

    message("Completed population analysis ", Sys.time())

    createSummaryExcelFromCumulativeBedFile(resultFolder = resultFolder, probeFeaturesOriginal = probeFeaturesOriginal, sampleSheet = sampleSheet)
    # xPlot <- colnames(dataToGraph)[-c('SAMPLENAME')] plot( xPlot, dataToGraph[2,2:dim(dataToGraph)[2]] ,type = 'l',col = 'red') for (z in 3:dim(dataToGraph)[1]) { lines(xPlot, dataToGraph[z,2:dim(dataToGraph)[2]], col = 'green') }
    end_time <- Sys.time()
    time_taken <- (end_time - start_time)
    message("Completed population with Excel summary",Sys.time(), " Time taken: ", time_taken)
}
