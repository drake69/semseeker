#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sample_sheet dataframe with at least a column Sample_ID to identify samples
#' @param methylation_data matrix of methylation data
#' @param bonferroni_threshold = 0.05 #threshold to define which pValue
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param iqrTimes how many times below the first quartile and over the third quartile the interqauartile is "added" to define the outlier
#' @param resultFolder wherte the result will be saved
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#'
semseeker <- function(sample_sheet,
                      methylation_data,
                      resultFolder,
                      bonferroni_threshold = 0.05,
                      maxResources = 90,
                      iqrTimes = 3 ) {

  envir <- init_env(resultFolder, maxResources)


  # set digits to 22
  withr::local_options(list(digits = 22))
  sliding_window_size <- 11


  methylation_data <- as.data.frame(methylation_data)
  methDataTemp <- data.frame( PROBE= rownames(methylation_data), methylation_data)
  methDataTemp <- methDataTemp[with(methDataTemp, order(methDataTemp$PROBE)), ]
  methylation_data <- methDataTemp[, -c(1)]

  rm(methDataTemp)
  message("INFO: I will work on:", nrow(methylation_data), " PROBES.")

  if(nrow(methylation_data) == 485512)
    message("INFO:seems a 450k dataset.")

  if(nrow(methylation_data) == 27578)
    message("INFO:seems a 27k dataset.")

  if(nrow(methylation_data) == 866562)
    message("INFO:seems an EPIC dataset.")

  PROBES <- PROBES[(PROBES$PROBE %in% rownames(methylation_data)),]
  methylation_data <- methylation_data[rownames(methylation_data) %in% PROBES$PROBE, ]
  methylation_data <- methylation_data[ order(rownames(methylation_data)), ]

  #PROBES <- sortByCHRandSTART(PROBES)
  if (!test_match_order(row.names(methylation_data), PROBES$PROBE)) {
    stop("Wrong order matching Probes and Methylation data!", Sys.time())
  }

  # if(!is.null(inferenceDetails))
  # {
  #   covariates <- unlist(t(strsplit( gsub(" ","",inferenceDetails$covariates),split = "+", fixed = T)))
  #
  #   if ( length(covariates)>0 & !( covariates  %in% colnames(sample_sheet)))
  #   {
  #     message(covariates[!(covariates%in% colnames(sample_sheet))])
  #     stop("The covariates value are not available in the sample sheet!")
  #   }
  # }


  populationCheckResult <- populationCheck(sample_sheet, methylation_data)
  if(!is.null(populationCheckResult))
  {
    stop(populationCheckResult)
  }

  # reference population
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(methylation_data), methylation_data[, referencePopulationSampleSheet$Sample_ID])

  #
  # methylation_data <- data.frame(PROBE = row.names(methylation_data), methylation_data[ , which(!(colnames(methylation_data)%in%referencePopulationSampleSheet$Sample_ID))]  )


  if (plyr::empty(referencePopulationMatrix) |
    dim(referencePopulationMatrix)[2] < 2) {
    message("Empty methylation_data ", Sys.time())
    stop("Empty methylation_data ")
  }

  populationControlRangeBetaValues <- rangeBetaValuePerProbeAsNormalizedDistribution(referencePopulationMatrix, iqrTimes)

  utils::write.table(x = populationControlRangeBetaValues, file = file_path_build(envir$resultFolderData ,"beta_thresholds","csv"), sep=";")

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)

  populations <-  c("Reference","Control","Case")
  for (populationName in populations) {

    #
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == populationName, ]
    populationMatrixColumns <- colnames(methylation_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: Population ",populationName, " is empty, probably the samples of this group are present in another group ? ", Sys.time())
      next
    }

    # browser()
    resultPopulation <- analize_population(
      envir = envir,
      methylation_data = methylation_data[, populationMatrixColumns] ,
      sliding_window_size = sliding_window_size,
      beta_superior_thresholds = populationControlRangeBetaValues$beta_superior_thresholds,
      beta_inferior_thresholds = populationControlRangeBetaValues$beta_inferior_thresholds,
      sample_sheet = populationSampleSheet,
      beta_medians = populationControlRangeBetaValues$betaMedianValues,
      bonferroni_threshold = bonferroni_threshold,
      probe_features = PROBES
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
  samplesID <- sample_sheet$Sample_ID
  sample_sheet <- sample_sheet[, !(colnames(sample_sheet) %in% colnames(resultSampleSheet))]
  sample_sheet$Sample_ID <- samplesID

  sample_sheet <- merge(sample_sheet, resultSampleSheet, by.x="Sample_ID", by.y="Sample_ID", all.x=TRUE)
  rm(methylation_data)
  gc()
  utils::write.csv2(sample_sheet, file.path(envir$resultFolderData , "sample_sheet_result.csv"))

  if(length(sample_sheet$Sample_Group=="Reference")>0){
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

