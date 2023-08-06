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

  probe_features <- probe_features_get("PROBE")
  message("DEBUG: ", Sys.time(), " loaded probe_features: PROBES")
  probe_features <- probe_features[(probe_features$PROBE %in% rownames(methylation_data)),]
  methylation_data <- methylation_data[rownames(methylation_data) %in% probe_features$PROBE, ]
  methylation_data <- methylation_data[ order(rownames(methylation_data)), ]

  # probe_features <- sort_by_chr_and_start(probe_features)
  if (!test_match_order(row.names(methylation_data), probe_features$PROBE)) {
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


  sample_group_checkResult <- sample_group_check(sample_sheet, methylation_data)
  if(!is.null(sample_group_checkResult))
  {
    stop(sample_group_checkResult)
  }

  beta_save(methylation_data, sample_sheet, batch_id)

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
  variables_to_export <- c( "ssEnv", "sample_sheet", "methylation_data", "analyze_population", "sliding_window_size",
    "populationControlRangeBetaValues", "bonferroni_threshold", "PROBES", "create_multiple_bed","probe_features")
  resultSampleSheet <- foreach::foreach(i = 1:length(ssEnv$keys_sample_groups[,1]), .combine = rbind, .export = variables_to_export ) %dorng%
  # for (i in 1:length(ssEnv$keys_sample_groups[,1]))
  {

    #
    sample_group <- ssEnv$keys_sample_groups[i,1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == sample_group, ]
    populationMatrixColumns <- colnames(methylation_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: ", Sys.time(), "  Population ",sample_group, " is empty, probably the samples of this group are present in another group ? ", Sys.time())
    }
    else
    {
      resultPopulation <- analyze_population(
        methylation_data = methylation_data[, populationMatrixColumns],
        sliding_window_size = sliding_window_size,
        sample_sheet = populationSampleSheet,
        beta_thresholds = populationControlRangeBetaValues,
        bonferroni_threshold = bonferroni_threshold,
        probe_features = probe_features
      )

      resultPopulation <- as.data.frame(resultPopulation)
      resultPopulation$Sample_Group <- sample_group
      create_multiple_bed(resultPopulation)

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
