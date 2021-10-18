#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sampleSheet dataframe with at least a column Sample_ID to identify
#' samples
#' @param methylationData matrix of methylation data
#' @param resultFolder folder to save computed epimutations and bedgraphs files.
#' @param bonferroniThreshold = 0.05 #threshold to define which pValue
#' @param covariates vector of column name from sample sheet to use as covariates
#' accept for lesions definition
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#'
semseeker <- function(sampleSheet,
                      methylationData,
                      resultFolder,
                      bonferroniThreshold = 0.05,
                      covariates = NULL) {

  sampleSheet <- as.data.frame(sampleSheet)
  # set digits to 22
  withr::local_options(list(digits = 22))
  slidingWindowSize <- 11

  if (resultFolder != "" && !dir.exists(resultFolder)) {
    dir.create(resultFolder)
  }

  methDataTemp <- data.frame( PROBE= rownames(methylationData), methylationData)
  methDataTemp <- methDataTemp[with(methDataTemp, order(PROBE)), ]

  methylationData <- methDataTemp[, -c(1)]

  rm(methDataTemp)
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

  logFolder <- paste0(tempdir(),"/log")

  populationControlRangeBetaValues <- NULL

  needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sampleSheet))]

  if (length(missedColumns) > 0) {
    stop(" Lost following columns ", missedColumns," ",Sys.time(), "Especting a column with name Sample_Group and possible values Reference,  Control and Case")
  }

  if (logFolder != "" && !dir.exists(logFolder)) {
    dir.create(logFolder)
  }

  populationGroups <- as.factor(sampleSheet$Sample_Group)
  # expected R as reference, S as Study, C as Control

  # reference population
  populations <-  c("Reference","Control","Case")
  sampleSheet$Sample_Group <- R.utils::toCamelCase(tolower(sampleSheet$Sample_Group), capitalize=TRUE)
  sampleSheet$Sample_Group <- as.factor(sampleSheet$Sample_Group)
  matchedPopulation <- levels(sampleSheet$Sample_Group) %in% populations
  if (is.element(FALSE, matchedPopulation)) {
    stop("File:",sampleSheetPath, " Sample_Group should contain only: Reference, Control, Case")
  }


  # reference population
  referencePopulationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(methylationData), methylationData[, referencePopulationSampleSheet$Sample_ID])

  # browser()
  # methylationData <- data.frame(PROBE = row.names(methylationData), methylationData[ , which(!(colnames(methylationData)%in%referencePopulationSampleSheet$Sample_ID))]  )


  if (plyr::empty(referencePopulationMatrix) |
    dim(referencePopulationMatrix)[2] < 2) {
    message("Empty methylationData ", Sys.time())
    stop("Empty methylationData ")
  }

  populationControlRangeBetaValues <- rangeBetaValuePerProbeAsNormalizedDistribution(referencePopulationMatrix, iqrTimes = 3)
  write.table(x = populationControlRangeBetaValues, file = file.path(resultFolder, "beta_thresholds.csv"), sep=";")

  for (i in 1:3) {

    # browser()
    # i <- 2
    populationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == populations[i], ]
    populationName <- populations[[i]]

    populationMatrixColumns <- colnames(methylationData[, populationSampleSheet$Sample_ID])
    # populationMatrixToAnalyze <- data.frame(PROBE = row.names(methylationData), methylationData[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: Population ",populationName, " is empty ", Sys.time())
      next
    }

    # extractEpiMutations(values = populationMatrixToAnalyze,resultFolder = resultFolder, thresholds = populationControlRangeBetaValues,
                        # populationName = populationSampleSheet$Sample_Group, probeFeatures= PROBES)

    analizePopulation(
      methylationData = methylationData,
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
    # rm(populationMatrixToAnalyze)
  }

  case_summary <- read.csv(file.path(resultFolder, "/Case/summary.csv"))
  ctrl_summary <- read.csv(file.path(resultFolder, "/Control/summary.csv"))
  reference_summary <-  read.csv(file.path(resultFolder, "/Reference/summary.csv"))
  studySummary <- rbind(case_summary, ctrl_summary, reference_summary)

  studySummary$MUTATIONS_BOTH <- studySummary$MUTATIONS_HYPER + studySummary$MUTATIONS_HYPO
  studySummary$LESIONS_BOTH <- studySummary$LESIONS_HYPER + studySummary$LESIONS_HYPO

  write.csv2(studySummary, file.path(resultFolder, "sample_sheet_result.csv"))


  populations <- c("Reference","Control","Case")
  mergeMultipleBed(
    populations,
    figures = c("METHYLATION"),
    anomalies = c("DELTAS"),
    fileExtension = ".bedgraph",
    resultFolder = resultFolder,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID", "VALUE")
  )

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")
  mergeMultipleBed(
    populations,
    figures,
    anomalies,
    fileExtension = ".bed",
    resultFolder = resultFolder,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID")
  )


  populations <- c("Reference","Control","Case")
  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")

  subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probesPrefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"
  try(
    createExcelPivot (logFolder, resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  )

  geneBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnID = c(3) ,resultFolder)
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE", groupColumnID = c(1) ,resultFolder)
  try(
    createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_PARTS", groupColumnID = c(1,3) ,resultFolder)
  )


  probesPrefix <- "PROBES_Island_"
  subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  mainGroupLabel <- "ISLAND"
  subGroupLabel <- "RELATION_TO_CPGISLAND"
  try(
    createExcelPivot (logFolder, resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  )

  islandBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "RELATION_TO_CPGISLAND", groupColumnID = 3 ,resultFolder)
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "ISLAND", groupColumnID = 1 ,resultFolder)
  try(
    createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "ISLAND_PARTS", groupColumnID = c(1,3) ,resultFolder)
  )

  subGroups <- c("DMR")
  probesPrefix = "PROBES_DMR_"
  mainGroupLabel =  "DMR"
  subGroupLabel="GROUP"
  try(
    createExcelPivot (logFolder, resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
  )

  dmrBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )
  createHeatmap(inputBedDataFrame =  dmrBed,anomalies = anomalies, groupLabel = mainGroupLabel, groupColumnID = 1 ,resultFolder)

  colnames(geneBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(dmrBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(islandBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")

  totalBed <- rbind(geneBed, dmrBed, islandBed, stringsAsFactors = TRUE)
  createHeatmap(inputBedDataFrame =  totalBed,anomalies = anomalies, groupLabel = "GENOMIC_AREA", groupColumnID = 3 ,resultFolder)

  rm(populationControlRangeBetaValues)

  message("Starting inference Analysis.")
  # inferenceAnalysis(studySummary = studySummary, resultFolder = resultFolder, logFolder= logFolder, family="gaussian", covariates= covariates, transformation = "log")
  fileName <- inferenceAnalysis(studySummary = studySummary, resultFolder = resultFolder, logFolder= logFolder, family="binomial", covariates= covariates, transformation = "log")
  geneontology_analysis_webgestalt(resultFolder = resultFolder, fileName = fileName)
  # inferenceAnalysis(studySummary = studySummary, resultFolder = resultFolder, logFolder= logFolder, family="poisson", covariates= covariates)
  # inferenceAnalysis(studySummary = studySummary, resultFolder = resultFolder, logFolder= logFolder, family="gaussian", covariates= covariates)
  # inferenceAnalysis(studySummary = studySummary, resultFolder = resultFolder, logFolder= logFolder, family="binomial", covariates= covariates)
  euristic_analysis_webgestalt(resultFolder = resultFolder)
  message("Job Completed !")
}

