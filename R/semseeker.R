#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sample_sheet dataframe with at least a column Sample_ID to identify samples
#' @param methylation_data matrix of methylation data
#' @param bonferroni_threshold = 0.05 #threshold to define which pValue adjusted to define an epilesion
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param iqrTimes how many times below the first quartile and over the third quartile the interqauartile is "added" to define the outlier
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param result_folder where the result will be saved
#'
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#'
semseeker <- function(sample_sheet,
                      methylation_data,
                      result_folder,
                      bonferroni_threshold = 0.05,
                      maxResources = 90,
                      iqrTimes = 3 ,
                      parallel_strategy ="multisession",
                      ... ) {

  envir <- init_env(result_folder, maxResources, parallel_strategy = parallel_strategy, ...)

  methylation_data <- stats::na.omit(methylation_data)
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

  #PROBES <- sort_by_chr_and_start(PROBES)
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


  population_checkResult <- population_check(sample_sheet, methylation_data, envir)
  if(!is.null(population_checkResult))
  {
    stop(population_checkResult)
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

  populationControlRangeBetaValues <- range_beta_values(referencePopulationMatrix, iqrTimes)

  utils::write.table(x = populationControlRangeBetaValues, file = file_path_build(envir$result_folderData ,"beta_thresholds","csv"), sep=";")

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)

  for (i in 1:length(envir$keys_populations[,1])) {

    #
    populationName <- envir$keys_populations[i,1]
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

    create_multiple_bed(envir, populationSampleSheet)

    resultPopulation <- as.data.frame(resultPopulation)
    # if(nrow(resultPopulation) != nrow(populationSampleSheet) )
    #   browser()

    if(!exists("resultSampleSheet"))
      resultSampleSheet <- resultPopulation
    else
      resultSampleSheet <- rbind(resultSampleSheet, resultPopulation)

    rm(populationSampleSheet)
  }

  sample_sheet <- as.data.frame(sample_sheet)
  resultSampleSheet <- as.data.frame(resultSampleSheet)
  samplesID <- sample_sheet$Sample_ID
  sample_sheet <- sample_sheet[, !(colnames(sample_sheet) %in% colnames(resultSampleSheet))]
  sample_sheet$Sample_ID <- samplesID

  sample_sheet <- merge(sample_sheet, resultSampleSheet, by.x="Sample_ID", by.y="Sample_ID", all.x=TRUE)
  rm(methylation_data)
  gc()
  utils::write.csv2(sample_sheet, file.path(envir$result_folderData , "sample_sheet_result.csv"))

  if(length(sample_sheet$Sample_Group=="Reference")>0)
    populations <- c("Reference","Control","Case")
  else
    populations <- c("Control","Case")

  figures <- envir$keys_figures[,1]
  anomalies <- envir$keys_anomalies[,1]

  if(sum(envir$keys_metaareas[,1]=="CHR")==1)
  {
    subGroups <- c("")
    probes_prefix = "PROBES"
    mainGroupLabel =  "CHR"
    subGroupLabel= "GROUP"

    create_excel_pivot (envir=envir, populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
    chrBed <- annotate_bed(envir=envir,populations ,figures ,anomalies ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( envir=envir,inputBedDataFrame =  chrBed,anomalies = anomalies, file_prefix = "CHR", groupColumnLabels = c("CHR"))
  }

  if(sum(envir$keys_metaareas[,1]=="GENE")==1)
  {
    subGroups <- envir$gene_subareas[,1]
    probes_prefix = "PROBES_Gene_"
    mainGroupLabel =  "GENE"
    subGroupLabel="GROUP"

    create_excel_pivot (envir=envir, populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

    geneBed <- annotate_bed(envir=envir,populations ,figures ,anomalies ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, file_prefix = "GENE_AREA", groupColumnLabels = c("GROUP"))
    create_heatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, file_prefix = "GENE", groupColumnLabels = c("GENE"))
    create_heatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, file_prefix = "GENE_PARTS", groupColumnLabels = c("GENE","GROUP"))
  }


  if(sum(envir$keys_metaareas[,1]=="ISLAND")==1)
  {
    probes_prefix <- "PROBES_Island_"
    subGroups <- envir$island_subareas[,1]
    mainGroupLabel <- "ISLAND"
    subGroupLabel <- "RELATION_TO_CPGISLAND"
    create_excel_pivot (envir=envir, populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)

    islandBed <- annotate_bed(envir=envir,populations ,figures ,anomalies ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( envir=envir,inputBedDataFrame =  islandBed,anomalies = anomalies, file_prefix = "RELATION_TO_CPGISLAND", groupColumnLabels = "RELATION_TO_CPGISLAND")
    create_heatmap( envir=envir,inputBedDataFrame =  islandBed,anomalies = anomalies, file_prefix = "ISLAND", groupColumnLabels = "ISLAND")
    create_heatmap( envir=envir,inputBedDataFrame =  geneBed,anomalies = anomalies, file_prefix = "ISLAND_PARTS", groupColumnLabels = c("ISLAND","RELATION_TO_CPGISLAND"))
  }

  if(sum(envir$keys_metaareas[,1]=="DMR")==1)
  {
    subGroups <- c("DMR")
    probes_prefix = "PROBES_DMR_"
    mainGroupLabel =  "DMR"
    subGroupLabel="GROUP"
    create_excel_pivot (envir=envir,populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)

    dmrBed <- annotate_bed(envir=envir,populations ,figures ,anomalies ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( envir=envir,inputBedDataFrame =  dmrBed,anomalies = anomalies, file_prefix = mainGroupLabel, groupColumnLabels = c("DMR") )
  }

  # if (!is.null(geneBed))
  # {
  #    geneBed <- geneBed[,c("MAINGROUP","SAMPLEID","SUBGROUP","VALUE","FIGURE","ANOMALY","POPULATION")]
  # }
  #
  # if (!is.null(dmrBed))
  # {
  #   dmrBed <- dmrBed[,c("MAINGROUP","SAMPLEID","SUBGROUP","VALUE","FIGURE","ANOMALY","POPULATION")]
  # }
  #
  # if (!is.null(islandBed))
  # {
  #   islandBed <- islandBed[,c("MAINGROUP","SAMPLEID","SUBGROUP","VALUE","FIGURE","ANOMALY","POPULATION")]
  # }
  #
  # totalBed <- rbind(geneBed, dmrBed, islandBed)
  # if (!is.null(totalBed) && nrow(totalBed)>0)
  #   create_heatmap( envir=envir,inputBedDataFrame =  totalBed,anomalies = anomalies, file_prefix = "GENOMIC_AREA", groupColumnLabels = 3)
  #
  # rm(populationControlRangeBetaValues)

  # message("Starting inference Analysis.")
  # inferenceAnalysis(envir$result_folderData = envir$result_folderData, envir$logFolder= envir$logFolder, inferenceDetails)
  # future::autoStopCluster(computationCluster)
  # doFuture::stopImplicitCluster()
  gc()
  # geneontology_analysis_webgestalt(envir$result_folderData = envir$result_folderData, fileName = fileName)
  # euristic_analysis_webgestalt(envir$result_folderData = envir$result_folderData)
  message("Job Completed !")

  future::plan( future::sequential)
}

