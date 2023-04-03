analyze_batch <- function(envir, methylation_data, sample_sheet, sliding_window_size, bonferroni_threshold,iqrTimes, batch_id)
{

  # browser()
  methylation_data <- as.data.frame(methylation_data)
  methDataTemp <- data.frame( "PROBE"= rownames(methylation_data), methylation_data)
  methDataTemp <- methDataTemp[with(methDataTemp, order(methDataTemp$PROBE)), ]
  methylation_data <- methDataTemp[, -c(1)]

  rm(methDataTemp)
  message("INFO: ", Sys.time(), " working on batch:", batch_id)
  message("INFO: ", Sys.time(), " I will work on:", nrow(methylation_data), " PROBES.")

  if(nrow(methylation_data) == 485512)
    message("INFO: ", Sys.time(), " seems a 450k dataset.")

  if(nrow(methylation_data) == 27578)
    message("INFO: ", Sys.time(), " seems a 27k dataset.")

  if(nrow(methylation_data) == 866562)
    message("INFO: ", Sys.time(), " seems an EPIC dataset.")

  probes <- semseeker::PROBES_CHR_CHR
  message("DEBUG: ", Sys.time(), " loaded probes: PROBES_CHR_CHR")
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


  if (plyr::empty(referencePopulationMatrix) ||
      ncol(referencePopulationMatrix) < 2) {
    message("ERROR: ", Sys.time(), " Empty methylation_data ", Sys.time())
    stop("INFO: ", Sys.time(), " Empty methylation_data ")
  }

  populationControlRangeBetaValues <- as.data.frame(range_beta_values(referencePopulationMatrix, iqrTimes))

  # utils::write.table(x = populationControlRangeBetaValues, file = file_path_build(envir$result_folderData ,c(batch_id, "beta_thresholds","csv")), sep=";")
  fst::write.fst(x = populationControlRangeBetaValues, path = file_path_build(envir$result_folderData ,c(batch_id, "beta_thresholds"),"fst"))

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)

  i <- 0
  # variables_to_export <- c( "envir", "sample_sheet", "methylation_data", "analize_population", "sliding_window_size", "populationControlRangeBetaValues", "bonferroni_threshold", "PROBES", "create_multiple_bed")
  # resultSampleSheet <- foreach::foreach(i = 1:length(envir$keys_populations[,1]), .combine = rbind, .export = variables_to_export ) %dorng%
  for (i in 1:length(envir$keys_populations[,1]))
  {

    #
    populationName <- envir$keys_populations[i,1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == populationName, ]
    populationMatrixColumns <- colnames(methylation_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      message("WARNING: ", Sys.time(), "  Population ",populationName, " is empty, probably the samples of this group are present in another group ? ", Sys.time())
    }
    else
    {
      # browser()
      resultPopulation <- analize_population(
        envir = envir,
        methylation_data = methylation_data[, populationMatrixColumns] ,
        sliding_window_size = sliding_window_size,
        beta_superior_thresholds = populationControlRangeBetaValues$beta_superior_thresholds,
        beta_inferior_thresholds = populationControlRangeBetaValues$beta_inferior_thresholds,
        sample_sheet = populationSampleSheet,
        beta_medians = populationControlRangeBetaValues$beta_median_values,
        bonferroni_threshold = bonferroni_threshold,
        probe_features = probes
      )

      resultPopulation <- as.data.frame(resultPopulation)
      resultPopulation$Sample_Group <- populationName
      create_multiple_bed(envir, resultPopulation)

      # resultPopulation
      # resultPopulation
      # # if(nrow(resultPopulation) != nrow(populationSampleSheet) )
      # #   browser()

      if(!exists("resultSampleSheet"))
        resultSampleSheet <- resultPopulation
      else
        resultSampleSheet <- plyr::rbind.fill(resultSampleSheet, resultPopulation)

      # rm(populationSampleSheet)
      gc()

    }
  }

  resultSampleSheet <- create_deltaq(envir, resultSampleSheet)

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
