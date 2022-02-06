#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sampleSheet dataframe with at least a column Sample_ID to identify
#' samples
#' @param methylationData matrix of methylation data
#' @param bonferroniThreshold = 0.05 #threshold to define which pValue
#' @param inferenceDetails dataframe of details to calculate inferential statistics accept for lesions definition
#' @param maxResources percentage of how many available cores will be used default 90%, rounded to the lowest integer
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#'
semseeker <- function(sampleSheet,
                      methylationData,
                      resultFolder,
                      bonferroniThreshold = 0.05,
                      inferenceDetails = NULL,
                      maxResources = 90) {


  init_env(resultFolder, maxResources)


  message(folderLog)
  sampleSheet <- as.data.frame(sampleSheet)
  sampleSheet <- sampleSheet[,!(colnames(sampleSheet) %in% c("Probes_Count", "MUTATIONS_HYPER", "LESIONS_HYPER", "MUTATIONS_HYPO", "LESIONS_HYPO", "MUTATIONS_BOTH", "LESIONS_BOTH"))]
  methylationData <- as.data.frame(methylationData)
  # set digits to 22
  withr::local_options(list(digits = 22))
  slidingWindowSize <- 11

  methDataTemp <- data.frame( PROBE= rownames(methylationData), methylationData)
  methDataTemp <- methDataTemp[with(methDataTemp, order(PROBE)), ]

  methylationData <- methDataTemp[, -c(1)]

  rm(methDataTemp)
  if(dim(methylationData)[1] <= 485512)
    message("WARNING: Actually semseeker works with EPIC, with 450K and 27K some probes will be lost.")

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

  populationControlRangeBetaValues <- NULL

  needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sampleSheet))]

  if (length(missedColumns) > 0) {
    stop(" Lost following columns ", missedColumns," ",Sys.time(), "Especting a column with name Sample_Group and possible values Reference,  Control and Case")
  }

  if (sum(colnames(methylationData) %in% sampleSheet$Sample_ID)!=ncol(methylationData))
  {
    stop("The methylation data has not the column's names as expected from the sample sheet in column Sample_ID!")
  }

  if(!is.null(inferenceDetails))
  {
    covariates <- unlist(t(strsplit( gsub(" ","",inferenceDetails$covariates),split = "+", fixed = T)))

    if ( length(covariates)>0 & !( covariates  %in% colnames(sampleSheet)))
    {
      message(covariates[!(covariates%in% colnames(sampleSheet))])
      stop("The covariates value are not available in the sample sheet!")
    }
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

  #
  # methylationData <- data.frame(PROBE = row.names(methylationData), methylationData[ , which(!(colnames(methylationData)%in%referencePopulationSampleSheet$Sample_ID))]  )


  if (plyr::empty(referencePopulationMatrix) |
    dim(referencePopulationMatrix)[2] < 2) {
    message("Empty methylationData ", Sys.time())
    stop("Empty methylationData ")
  }


  populationControlRangeBetaValues <- rangeBetaValuePerProbeAsNormalizedDistribution(referencePopulationMatrix, iqrTimes = 3)

  write.table(x = populationControlRangeBetaValues, file = file_path_build(resultFolderData ,"beta_thresholds","csv"), sep=";")

  # remove duplicated samples due to the reference population
  referenceSamples <- sampleSheet[sampleSheet$Sample_Group == "Reference",]
  otherSamples <- sampleSheet[sampleSheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sampleSheet <- rbind(otherSamples, referenceSamples)

  for (populationName in populations) {

    #
    populationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == populationName, ]
    populationMatrixColumns <- colnames(methylationData[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: Population ",populationName, " is empty ", Sys.time())
      next
    }

    # browser()
    resultPopulation <- analizePopulation(
      methylationData = methylationData[, populationMatrixColumns] ,
      slidingWindowSize = slidingWindowSize,
      betaSuperiorThresholds = populationControlRangeBetaValues$betaSuperiorThresholds,
      betaInferiorThresholds = populationControlRangeBetaValues$betaInferiorThresholds,
      sampleSheet = populationSampleSheet,
      betaMedians = populationControlRangeBetaValues$betaMedianValues,
      bonferroniThreshold = bonferroniThreshold,
      probeFeatures = PROBES
    )

    createMultipleBed(populationSampleSheet)

    resultPopulation <- as.data.frame(resultPopulation)
    # if(nrow(resultPopulation) != nrow(populationSampleSheet) )
    #   browser()

    if(!exists("resultSampleSheet"))
      resultSampleSheet <- resultPopulation
    else
      resultSampleSheet <- rbind(resultSampleSheet, resultPopulation)

    rm(populationSampleSheet)
  }

  # browser()
  samplesID <- sampleSheet$Sample_ID
  sampleSheet <- sampleSheet[, !(colnames(sampleSheet) %in% colnames(resultSampleSheet))]
  sampleSheet$Sample_ID <- samplesID

  sampleSheet <- merge(sampleSheet, resultSampleSheet, by.x="Sample_ID", by.y="Sample_ID", all.x=TRUE)
  rm(methylationData)
  gc()
  write.csv2(sampleSheet, file.path(resultFolderData , "sample_sheet_result.csv"))

  if(length(sampleSheet$Sample_Group=="Reference")>0){
    populations <- c("Reference","Control","Case")
  }  else
    populations <- c("Control","Case")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")

  subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probesPrefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"

  createExcelPivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probesPrefix =   probesPrefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  geneBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnID = c(3))
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE", groupColumnID = c(1) )
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_PARTS", groupColumnID = c(1,3))


  probesPrefix <- "PROBES_Island_"
  subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  mainGroupLabel <- "ISLAND"
  subGroupLabel <- "RELATION_TO_CPGISLAND"
  createExcelPivot ( populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)

  islandBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "RELATION_TO_CPGISLAND", groupColumnID = 3)
  createHeatmap(inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "ISLAND", groupColumnID = 1)
  createHeatmap(inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "ISLAND_PARTS", groupColumnID = c(1,3))

  subGroups <- c("DMR")
  probesPrefix = "PROBES_DMR_"
  mainGroupLabel =  "DMR"
  subGroupLabel="GROUP"
  createExcelPivot (populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)

  dmrBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,)
  createHeatmap(inputBedDataFrame =  dmrBed,anomalies = anomalies, groupLabel = mainGroupLabel, groupColumnID = 1 )

  colnames(geneBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(dmrBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  colnames(islandBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")

  totalBed <- rbind(geneBed, dmrBed, islandBed, stringsAsFactors = TRUE)
  createHeatmap(inputBedDataFrame =  totalBed,anomalies = anomalies, groupLabel = "GENOMIC_AREA", groupColumnID = 3)

  rm(populationControlRangeBetaValues)

  # message("Starting inference Analysis.")
  # inferenceAnalysis(resultFolder = resultFolder, folderLog= folderLog, inferenceDetails)
  parallel::stopCluster(computationCluster)
  doParallel::stopImplicitCluster()
  gc()
  # geneontology_analysis_webgestalt(resultFolder = resultFolder, fileName = fileName)
  # euristic_analysis_webgestalt(resultFolder = resultFolder)
  message("Job Completed !")


}


semseeker_pheno <- function(sampleSheet,
                      methylationData,
                      resultFolder,
                      bonferroniThreshold = 0.05,
                      inferenceDetails = NULL,
                      pheno_term,
                      onlySeed = TRUE) {

  probesFilter <- probes_go_association_phenolizer(pheno_term, onlySeed = TRUE )

  # browser()
  methylationData <- methylationData[ rownames(methylationData) %in% probesFilter, ]

  semseeker(sampleSheet = sampleSheet,
            methylationData = methylationData,
            resultFolder = resultFolder,
            inferenceDetails = inferenceDetails)



}
