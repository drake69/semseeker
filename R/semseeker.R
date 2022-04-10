#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sampleSheet dataframe with at least a column Sample_ID to identify samples
#' @param methylationData matrix of methylation data
#' @param bonferroniThreshold = 0.05 #threshold to define which pValue
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param iqrTimes how many times below the first quartile and over the third quartile the interqauartile is "added" to define the outlier
#' @param resultFolder wherte the result will be saved
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#'
semseeker <- function(sampleSheet,
                      methylationData,
                      resultFolder,
                      bonferroniThreshold = 0.05,
                      maxResources = 90,
                      iqrTimes = 3 ) {

  envir <- init_env(resultFolder, maxResources)


  # set digits to 22
  withr::local_options(list(digits = 22))
  slidingWindowSize <- 11


  methylationData <- as.data.frame(methylationData)
  methDataTemp <- data.frame( PROBE= rownames(methylationData), methylationData)
  methDataTemp <- methDataTemp[with(methDataTemp, order(methDataTemp$PROBE)), ]
  methylationData <- methDataTemp[, -c(1)]

  rm(methDataTemp)
  message("INFO: I will work on:", nrow(methylationData), " PROBES.")

  if(nrow(methylationData) == 485512)
    message("INFO:seems a 450k dataset.")

  if(nrow(methylationData) == 27578)
    message("INFO:seems a 27k dataset.")

  if(nrow(methylationData) == 866562)
    message("INFO:seems an EPIC dataset.")

  PROBES <- PROBES[(PROBES$PROBE %in% rownames(methylationData)),]
  methylationData <- methylationData[rownames(methylationData) %in% PROBES$PROBE, ]
  methylationData <- methylationData[ order(rownames(methylationData)), ]

  #PROBES <- sortByCHRandSTART(PROBES)
  if (!test_match_order(row.names(methylationData), PROBES$PROBE)) {
    stop("Wrong order matching Probes and Methylation data!", Sys.time())
  }

  # if(!is.null(inferenceDetails))
  # {
  #   covariates <- unlist(t(strsplit( gsub(" ","",inferenceDetails$covariates),split = "+", fixed = T)))
  #
  #   if ( length(covariates)>0 & !( covariates  %in% colnames(sampleSheet)))
  #   {
  #     message(covariates[!(covariates%in% colnames(sampleSheet))])
  #     stop("The covariates value are not available in the sample sheet!")
  #   }
  # }


  populationCheckResult <- populationCheck(sampleSheet, methylationData)
  if(!is.null(populationCheckResult))
  {
    stop(populationCheckResult)
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

  populationControlRangeBetaValues <- rangeBetaValuePerProbeAsNormalizedDistribution(referencePopulationMatrix, iqrTimes)

  utils::write.table(x = populationControlRangeBetaValues, file = file_path_build(envir$resultFolderData ,"beta_thresholds","csv"), sep=";")

  # remove duplicated samples due to the reference population
  referenceSamples <- sampleSheet[sampleSheet$Sample_Group == "Reference",]
  otherSamples <- sampleSheet[sampleSheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sampleSheet <- rbind(otherSamples, referenceSamples)

  populations <-  c("Reference","Control","Case")
  for (populationName in populations) {

    #
    populationSampleSheet <- sampleSheet[sampleSheet$Sample_Group == populationName, ]
    populationMatrixColumns <- colnames(methylationData[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: Population ",populationName, " is empty, probably the samples of this group are present in another group ? ", Sys.time())
      next
    }

    # browser()
    resultPopulation <- analizePopulation(
      envir = envir,
      methylationData = methylationData[, populationMatrixColumns] ,
      slidingWindowSize = slidingWindowSize,
      betaSuperiorThresholds = populationControlRangeBetaValues$betaSuperiorThresholds,
      betaInferiorThresholds = populationControlRangeBetaValues$betaInferiorThresholds,
      sampleSheet = populationSampleSheet,
      betaMedians = populationControlRangeBetaValues$betaMedianValues,
      bonferroniThreshold = bonferroniThreshold,
      probeFeatures = PROBES
    )

    createMultipleBed(envir, populationSampleSheet)

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
  utils::write.csv2(sampleSheet, file.path(envir$resultFolderData , "sample_sheet_result.csv"))

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

  createExcelPivot (envir=envir, populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probesPrefix =   probesPrefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  geneBed <- annotateBed(envir=envir,populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)
  createHeatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnIDs = c(3))
  createHeatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE", groupColumnIDs = c(1) )
  createHeatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "GENE_PARTS", groupColumnIDs = c(1,3))


  probesPrefix <- "PROBES_Island_"
  subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
  mainGroupLabel <- "ISLAND"
  subGroupLabel <- "RELATION_TO_CPGISLAND"
  createExcelPivot (envir=envir, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)

  islandBed <- annotateBed(envir=envir,populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)
  createHeatmap( envir=envir,inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "RELATION_TO_CPGISLAND", groupColumnIDs = 3)
  createHeatmap( envir=envir,inputBedDataFrame =  islandBed,anomalies = anomalies, groupLabel = "ISLAND", groupColumnIDs = 1)
  createHeatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, groupLabel = "ISLAND_PARTS", groupColumnIDs = c(1,3))

  subGroups <- c("DMR")
  probesPrefix = "PROBES_DMR_"
  mainGroupLabel =  "DMR"
  subGroupLabel="GROUP"
  createExcelPivot (envir=envir,populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)

  dmrBed <- annotateBed(envir=envir,populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)
  createHeatmap( envir=envir,inputBedDataFrame =  dmrBed,anomalies = anomalies, groupLabel = mainGroupLabel, groupColumnIDs = 1 )

  if (!is.null(geneBed))
  {
    colnames(geneBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  }

  if (!is.null(dmrBed))
  {
    colnames(dmrBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  }

  if (!is.null(islandBed))
  {
    colnames(islandBed) <- c("MAINGROUP","SAMPLEID","SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  }

  totalBed <- rbind(geneBed, dmrBed, islandBed)
  if (!is.null(totalBed) && nrow(totalBed)>0)
    createHeatmap( envir=envir,inputBedDataFrame =  totalBed,anomalies = anomalies, groupLabel = "GENOMIC_AREA", groupColumnIDs = 3)

  rm(populationControlRangeBetaValues)

  # message("Starting inference Analysis.")
  # inferenceAnalysis(envir$resultFolderData = envir$resultFolderData, envir$logFolder= envir$logFolder, inferenceDetails)
  # future::autoStopCluster(computationCluster)
  # doFuture::stopImplicitCluster()
  gc()
  # geneontology_analysis_webgestalt(envir$resultFolderData = envir$resultFolderData, fileName = fileName)
  # euristic_analysis_webgestalt(envir$resultFolderData = envir$resultFolderData)
  message("Job Completed !")


}

