#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sampleSheetPath dataframe with at least a column Sample_ID to identify
#' samples
#' @param methylationData matrix of methylation data
#' @param resultFolder folder to save computed epimutations and bedgraphs files.
#' @param bonferroniThreshold = 0.05 #threshold to define which pValue
#' accept for lesions definition
#' @return files into the result folder with pivot table and bedgraph.
#' @export
calculate <- function(sampleSheetPath,
                      methylationData,
                      resultFolder,
                      bonferroniThreshold = 0.05) {

  # set digits to 22
  withr::local_options(list(digits = 22))
  slidingWindowSize <- 11

  if (resultFolder != "" && !dir.exists(resultFolder)) {
    dir.create(resultFolder)
  }

  methDataTemp <- data.frame( PROBE= rownames(methylationData), methylationData)
  methDataTemp <- methDataTemp[with(methDataTemp, order(PROBE)), ]

  methylationData <- methDataTemp[, -c(1)]

  # if(dim(methylationData)[1] < 485512)
  #   PROBES <- PROBES[!is.na(PROBES$METH_450K),]
  #
  # if(dim(methylationData)[1] < 27578)
  #   PROBES <- PROBES[!is.na(PROBES$METH_27K),]

  PROBES <- PROBES[(PROBES$PROBE %in% rownames(methylationData)),]

  methylationData <- methylationData[rownames(methylationData) %in% PROBES$PROBE, ]

  #PROBES <- sortByCHRandSTART(PROBES)
  if (!test_match_order(row.names(methylationData), PROBES$PROBE)) {
    stop("Wrong order matching Probes and Methylation data!", Sys.time())
  }


  sampleSheet <- utils::read.csv(sampleSheetPath)
  logFolder <- paste0(tempdir(),"/log")

  populationControlRangeBetaValues <- NULL

  needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sampleSheet))]

  if (length(missedColumns) > 0) {
    stop("File:",sampleSheet, " Lost following columns ", missedColumns," ",Sys.time(), "Especting a column with name Sample_Group and possible values 0 as reference, 1 as Control 2 as Case")
  }

  if (logFolder != "" && !dir.exists(logFolder)) {
    dir.create(logFolder)
  }

  populationGroups <- as.factor(sampleSheet$Sample_Group)
  # expected R as reference, S as Study, C as Control

  # reference population
  populations <-  c("Reference","Control","Case")

  # reference population
  referencePopulationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(methylationData), methylationData[, referencePopulationSampleSheet$Sample_ID])

  if (plyr::empty(referencePopulationMatrix) |
    dim(referencePopulationMatrix)[2] < 2) {
    message("Empty methylationData ", Sys.time())
    stop("Empty methylationData ")
  }

  populationControlRangeBetaValues <- rangeBetaValuePerProbeAsNormalizedDistribution(referencePopulationMatrix, iqrTimes = 3)


  for (i in 1:3) {

    populationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == populations[i], ]
    populationName <- populations[[i]]

    populationMatrixToAnalyze <- data.frame(PROBE = row.names(methylationData), methylationData[, populationSampleSheet$Sample_ID])

    if (plyr::empty(populationMatrixToAnalyze) |
      dim(populationMatrixToAnalyze)[2] < 2) {
      message("WARNING: Population ",populationName, " is empty ", Sys.time())
      next
    }

    analizePopulation(
      populationMatrix = populationMatrixToAnalyze,
      slidingWindowSize = slidingWindowSize,
      resultFolder = resultFolder,
      logFolder = logFolder,
      betaSuperiorThresholds = populationControlRangeBetaValues$betaSuperiorThresholds,
      betaInferiorThresholds = populationControlRangeBetaValues$betaInferiorThresholds,
      sampleSheet = populationSampleSheet,
      betaMedians = populationControlRangeBetaValues$betaMedianValues,
      populationName = populationName,
      bonferroniThreshold = bonferroniThreshold,
      probeFeatures = PROBES
    )

    rm(populationSampleSheet)
    rm(populationMatrixToAnalyze)
  }


  populations <- c("Reference","Control","Case")
  mergeMultipleBed(
    populations,
    figures = c("METHYLATION"),
    anomalies = c("DELTAS"),
    fileExtension = ".bedgraph",
    resultFolder = resultFolder,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME", "VALUE")
  )

  figures <- c("HYPO", "HYPER")
  anomalies <- c("MUTATIONS","LESIONS")
  mergeMultipleBed(
    populations,
    figures,
    anomalies,
    fileExtension = ".bed",
    resultFolder = resultFolder,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  # createSummaryExcelFromCumulativeBedFile(resultFolder = resultFolder, probeFeatures = probeFeatures,
                                          # sampleSheet = sampleSheet)

  populations <- c("Reference","Control","Case")
  figures <- c("HYPO", "HYPER")
  anomalies <- c("MUTATIONS","LESIONS")

  subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd")
  probesPrefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"
  # regionDifferentialAnalysisPerGenomicArea(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  createChartFromMultipleBedGenericPerRegion(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel )
  # regionDifferentialAnalysis(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  geneBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnID = c(3) )
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE", groupColumnID = c(1) )
  try(
    createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_PARTS", groupColumnID = c(1,3) )
  )

  probesPrefix <- "PROBES_Island_"
  subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island")
  mainGroupLabel <- "ISLAND"
  subGroupLabel <- "RELATION_TO_CPGISLAND"
  islandBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  # regionDifferentialAnalysisPerGenomicArea(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  createChartFromMultipleBedGenericPerRegion(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel )
  # regionDifferentialAnalysis(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "RELATION_TO_CPGISLAND", groupColumnID = 3 )
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "ISLAND", groupColumnID = 1 )
  try(
    createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "ISLAND_PARTS", groupColumnID = c(1,3) )
  )

  subGroups <- c("DMR")
  probesPrefix = "PROBES_DMR_"
  mainGroupLabel =  "DMR"
  subGroupLabel="GROUP"
  dmrBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  # regionDifferentialAnalysisPerGenomicArea(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  createChartFromMultipleBedGenericPerRegion(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel )
  # regionDifferentialAnalysis(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  createHeatmap(inputBedDataFrame =  dmrBed,anomalies = anomalies, groupLabel = mainGroupLabel, groupColumnID = 1 )

  colnames(geneBed) <- c("MAINGROUP","SAMPLENAME","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(dmrBed) <- c("MAINGROUP","SAMPLENAME","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(islandBed) <- c("MAINGROUP","SAMPLENAME","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  totalBed <- rbind(geneBed, dmrBed, islandBed, stringsAsFactors = TRUE)
  createHeatmap(inputBedDataFrame =  totalBed,anomalies = anomalies, groupLabel = "GENOMIC_AREA", groupColumnID = 3 )

  rm(populationControlRangeBetaValues)
}
