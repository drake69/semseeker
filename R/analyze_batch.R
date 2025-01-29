analyze_batch <- function(signal_data, sample_sheet, batch_id)
{
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " working on batch:", batch_id)

  ssEnv <- get_session_info()

  # #
  signal_data <- as.data.frame(signal_data)
  get_meth_tech(signal_data)
  coverage_analysis(signal_data = signal_data)

  methDataTemp <- data.frame( "PROBE"= rownames(signal_data), signal_data)
  methDataTemp <- methDataTemp[with(methDataTemp, order(methDataTemp$PROBE)), ]
  signal_data <- methDataTemp[, -c(1)]
  rm(methDataTemp)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work on:", nrow(signal_data), " PROBES.")

  probe_features <- probe_features_get("PROBE")
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " loaded probe_features: PROBES")

  probe_features <- probe_features[(probe_features$PROBE %in% rownames(signal_data)),]
  signal_data <- signal_data[rownames(signal_data) %in% probe_features$PROBE, ]
  signal_data <- signal_data[ order(rownames(signal_data)), ]

  # probe_features <- sort_by_chr_and_start(probe_features)
  if (!test_match_order(row.names(signal_data), probe_features$PROBE)) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Wrong order matching Probes and Methylation data!")
    stop()
  }

  # if(!is.null(inferenceDetails))
  # {
  #   covariates <- unlist(t(strsplit( gsub(" ","",inferenceDetails$covariates),split = "+", fixed = T)))
  #
  #   if ( length(covariates)>0 & !( covariates  %in% colnames(sample_sheet)))
  #   {
  #     log_event(covariates[!(covariates%in% colnames(sample_sheet))])
  #     stop("The covariates value are not available in the sample sheet!")
  #   }
  # }


  sample_group_checkResult <- sample_group_check(sample_sheet, signal_data)
  if(!is.null(sample_group_checkResult))
  {
    stop(sample_group_checkResult)
  }

  signal_save(signal_data, sample_sheet, batch_id)

  # reference population
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(signal_data), signal_data[, referencePopulationSampleSheet$Sample_ID])

  #
  # signal_data <- data.frame(PROBE = row.names(signal_data), signal_data[ , which(!(colnames(signal_data)%in%referencePopulationSampleSheet$Sample_ID))]  )


  if (plyr::empty(referencePopulationMatrix) ||
      ncol(referencePopulationMatrix) < 2) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Empty signal_data ", format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  if(!ssEnv$signal_intrasample)
  {
    populationControlRangeBetaValues <- as.data.frame(signal_range_values(referencePopulationMatrix))

    # utils::write.table(x = gzfile(populationControlRangeBetaValues), file = file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds","csv"), add_gz = TRUE), sep=";")
    fst::write.fst(x = populationControlRangeBetaValues, path = file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"fst"), compress = 100)
  }
  else
  {
    populationControlRangeBetaValues <- NULL
  }

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)

  i <- 0
  variables_to_export <- c( "ssEnv", "sample_sheet", "signal_data", "analyze_population",
    "populationControlRangeBetaValues", "PROBES", "create_multiple_bed","probe_features")
  # resultSampleSheet <- foreach::foreach(i = 1:length(ssEnv$keys_sample_groups[,1]), .combine = rbind, .export = variables_to_export ) %dorng%
  for (i in 1:length(ssEnv$keys_sample_groups[,1]))
  {
    #
    sample_group <- ssEnv$keys_sample_groups[i,1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == sample_group, ]
    populationMatrixColumns <- colnames(signal_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), "  Population ",sample_group, " is empty, probably the samples of this group are present in another group ? ")
    }
    else
    {

      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "  Working on population ",sample_group, " with ", nrow(signal_data), " probes.")
      resultPopulation <- analyze_population(
        signal_data = signal_data[, populationMatrixColumns],
        sample_sheet = populationSampleSheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features = probe_features
      )

      resultPopulation <- as.data.frame(resultPopulation)
      resultPopulation$Sample_Group <- sample_group
      create_multiple_bed(resultPopulation)

      # resultPopulation
      # resultPopulation
      # # if(nrow(resultPopulation) != nrow(populationSampleSheet) )
      # #   #

      if(!exists("resultSampleSheet"))
        resultSampleSheet <- resultPopulation
      else
        resultSampleSheet <- plyr::rbind.fill(resultSampleSheet, resultPopulation)

      # rm(populationSampleSheet)

      # resultPopulation
    }
  }

  # resultSampleSheet <- as.data.frame(resultSampleSheet)
  resultSampleSheet <- deltaq_get(resultSampleSheet)
  resultSampleSheet <- deltarq_get(resultSampleSheet)

  resultSampleSheet <- deltap_get(resultSampleSheet)
  resultSampleSheet <- deltarp_get(resultSampleSheet)

  sample_sheet <- as.data.frame(sample_sheet)
  resultSampleSheet <- as.data.frame(resultSampleSheet)
  # fill na with zeros
  resultSampleSheet[is.na(resultSampleSheet)] <- 0

  samplesID <- sample_sheet$Sample_ID
  sample_sheet <- sample_sheet[, !(colnames(sample_sheet) %in% colnames(resultSampleSheet))]
  sample_sheet$Sample_ID <- samplesID
  sample_sheet$Batch_ID <- batch_id

  sample_sheet <- merge(sample_sheet, resultSampleSheet, by.x="Sample_ID", by.y="Sample_ID", all.x=TRUE)
  rm(signal_data)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "  Batch completed:", batch_id)
  return((sample_sheet))

}
