analyze_batch <- function(methylation_data, sample_sheet, sliding_window_size, bonferroni_threshold,iqrTimes, batch_id)
{
  message("INFO: ", Sys.time(), " working on batch:", batch_id)

  ssEnv <- get_session_info()

  # browser()
  methylation_data <- as.data.frame(methylation_data)
  get_meth_tech(methylation_data)
  coverage_analysis(methylation_data = methylation_data)

  methDataTemp <- data.frame( "PROBE"= rownames(methylation_data), methylation_data)
  methDataTemp <- methDataTemp[with(methDataTemp, order(methDataTemp$PROBE)), ]
  methylation_data <- methDataTemp[, -c(1)]
  rm(methDataTemp)

  message("INFO: ", Sys.time(), " I will work on:", nrow(methylation_data), " PROBES.")

  probes <- semseeker::PROBES
  message("DEBUG: ", Sys.time(), " loaded probes: PROBES")
  probes <- probes[(probes$PROBE %in% rownames(methylation_data)),]
  methylation_data <- methylation_data[rownames(methylation_data) %in% probes$PROBE, ]
  methylation_data <- methylation_data[ order(rownames(methylation_data)), ]

  # probes <- sort_by_chr_and_start(probes)
  if (!test_match_order(row.names(methylation_data), probes$PROBE)) {
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


  population_checkResult <- population_check(sample_sheet, methylation_data)
  if(!is.null(population_checkResult))
  {
    stop(population_checkResult)
  }

  # # save beta as rds and as pivot
  # sample_info <- sample_sheet[,c("Sample_ID","Sample_Group")]
  # colnames(sample_info) <- c("SAMPLEID","POPULATION")
  # methylation_data_to_save <- data.frame("SAMPLEID"=rownames(methylation_data), methylation_data[, sample_sheet$Sample_ID])
  # sample_info <- as.data.frame(t(sample_info))
  # colnames(sample_info) <- sample_info[1,]
  # sample_info <- cbind(data.frame("SAMPLEID"="POPULATION"), sample_info[-1,])
  # methylation_data_to_save <- rbind(sample_info, methylation_data_to_save )
  # pivot_subfolder <- dir_check_and_create(ssEnv$result_folderData, "BETA")
  # fileName <- file_path_build(pivot_subfolder,c("BETA","BETA","PROBE","PROBE"),"csv")
  # # fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
  # utils::write.table(methylation_data_to_save, fileName, row.names = T, col.names = T, sep=";")
  # rm(methylation_data_to_save)
  # rm(sample_info)
  # saveRDS(methylation_data,file_path_build(ssEnv$result_folderData, c(batch_id,"_methylation_data"),"rds"))

  # reference population
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(methylation_data), methylation_data[, referencePopulationSampleSheet$Sample_ID])

  #
  # methylation_data <- data.frame(PROBE = row.names(methylation_data), methylation_data[ , which(!(colnames(methylation_data)%in%referencePopulationSampleSheet$Sample_ID))]  )


  if (plyr::empty(referencePopulationMatrix) ||
      ncol(referencePopulationMatrix) < 2) {
    message("ERROR: ", Sys.time(), " Empty methylation_data ", Sys.time())
    stop("INFO: ", Sys.time(), " Empty methylation_data ")
  }

  populationControlRangeBetaValues <- as.data.frame(range_beta_values(referencePopulationMatrix, iqrTimes))

  # utils::write.table(x = populationControlRangeBetaValues, file = file_path_build(ssEnv$result_folderData ,c(batch_id, "beta_thresholds","csv")), sep=";")
  fst::write.fst(x = populationControlRangeBetaValues, path = file_path_build(ssEnv$result_folderData ,c(batch_id, "beta_thresholds"),"fst"))

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)

  i <- 0
  variables_to_export <- c( "ssEnv", "sample_sheet", "methylation_data", "analize_population", "sliding_window_size",
    "populationControlRangeBetaValues", "bonferroni_threshold", "PROBES", "create_multiple_bed","probes")
  resultSampleSheet <- foreach::foreach(i = 1:length(ssEnv$keys_populations[,1]), .combine = rbind, .export = variables_to_export ) %dorng%
  # for (i in 1:length(ssEnv$keys_populations[,1]))
  {

    #
    populationName <- ssEnv$keys_populations[i,1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == populationName, ]
    populationMatrixColumns <- colnames(methylation_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: ", Sys.time(), "  Population ",populationName, " is empty, probably the samples of this group are present in another group ? ", Sys.time())
    }
    else
    {
      resultPopulation <- analize_population(
        methylation_data = methylation_data[, populationMatrixColumns],
        sliding_window_size = sliding_window_size,
        sample_sheet = populationSampleSheet,
        beta_thresholds = populationControlRangeBetaValues,
        bonferroni_threshold = bonferroni_threshold,
        probe_features = probes
      )

      resultPopulation <- as.data.frame(resultPopulation)
      resultPopulation$Sample_Group <- populationName
      create_multiple_bed( resultPopulation)

      # resultPopulation
      # resultPopulation
      # # if(nrow(resultPopulation) != nrow(populationSampleSheet) )
      # #   browser()

      # if(!exists("resultSampleSheet"))
      #   resultSampleSheet <- resultPopulation
      # else
      #   resultSampleSheet <- plyr::rbind.fill(resultSampleSheet, resultPopulation)

      # rm(populationSampleSheet)
      gc()
      resultPopulation
    }
  }

  # resultSampleSheet <- as.data.frame(resultSampleSheet)
  resultSampleSheet <- deltaq_get(resultSampleSheet)
  resultSampleSheet <- deltarq_get(resultSampleSheet)

  sample_sheet <- as.data.frame(sample_sheet)
  resultSampleSheet <- as.data.frame(resultSampleSheet)
  samplesID <- sample_sheet$Sample_ID
  sample_sheet <- sample_sheet[, !(colnames(sample_sheet) %in% colnames(resultSampleSheet))]
  sample_sheet$Sample_ID <- samplesID
  sample_sheet$Batch_ID <- batch_id

  sample_sheet <- merge(sample_sheet, resultSampleSheet, by.x="Sample_ID", by.y="Sample_ID", all.x=TRUE)
  rm(methylation_data)

  message("INFO: ", Sys.time(), "  Batch completed:", batch_id)
  return((sample_sheet))

}
