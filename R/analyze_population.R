#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param signal_data whole matrix of data to analyze.
#' @param sample_sheet name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param probe_features probe_features detail from 27 to EPIC illumina dataset
#' @param signal_thresholds thresholds defined to calculate epimutations
#' @return files into the result folder with pivot table and bedgraph.
#' @importFrom doRNG %dorng%
#'
analyze_population <- function(signal_data, sample_sheet,signal_thresholds, probe_features) {

  ssEnv <- semseeker:::get_session_info()
  # #
  start_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " AnalyzePopulation warmingUP ")

  signal_data <- stats::na.omit(signal_data)

  ### get signal_values ########################################################
  sample_sheet <- sample_sheet[order(sample_sheet[, "Sample_ID"], decreasing = FALSE), ]
  existent_samples <- colnames(signal_data)
  sample_names <- sample_sheet$Sample_ID
  missed_samples <- setdiff(setdiff(sample_names, existent_samples), "PROBE")

  if (length(missed_samples) != 0) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " These samples data are missed: ", paste0(missed_samples, sep = " "))
  }

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " WarmedUP AnalyzePopulation")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Start population analysis")

  # progress_bar <- progress::progress_bar$new(
  #   format = paste("INFO: Performing population analysis [:bar] :percent eta: :eta"),
  #   total = nrow(sample_sheet),
  #   clear = FALSE,
  #   width= 60)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  variables_to_export <- c("sample_sheet", "signal_data", "semseeker:::analyze_single_sample", "ssEnv",
                            "signal_superior_thresholds","deltar_single_sample","signal_inferior_thresholds","iqr","signal_median_values",
                           "bt","bonferroni_threshold", "probe_features", "semseeker:::analyze_single_sample_both", "delta_single_sample", "progress_bar",
                           "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_single_sample",
                           "get_session_info")
  i <- 1

  if (!ssEnv$signal_intrasample)
  {
    signal_superior_thresholds <- signal_thresholds$signal_superior_thresholds
    signal_inferior_thresholds <- signal_thresholds$signal_inferior_thresholds
    iqr <- signal_thresholds$iqr
    signal_median_values <- signal_thresholds$signal_median_values
  }

  # for(i in 1:nrow(sample_sheet)) {
  summary_population <-  foreach::foreach(i =1:nrow(sample_sheet), .combine= plyr::rbind.fill , .export = variables_to_export) %dorng% {
    local_sample_detail <- sample_sheet[i,]

    ssEnv <- semseeker:::get_session_info()
    signal_values <- signal_data[, local_sample_detail$Sample_ID]
    if (ssEnv$signal_intrasample )
    {
      q <- stats::quantile(signal_values)
      q1 <- as.numeric(q[2])
      q3 <- as.numeric(q[4])
      y_med <- as.numeric(q[3])
      iqr <- stats::IQR(signal_values)
      iqrmult <- 3
      y_sup <- q3 + iqrmult * iqr
      y_inf <- q1 - iqrmult * iqr
      signal_superior_thresholds <- rep(y_sup,length(signal_values))
      signal_inferior_thresholds <- rep(y_inf,length(signal_values))
      signal_median_values <- rep(y_med,length(signal_values))
    }

    hyper_result <- semseeker:::analyze_single_sample( values = signal_values,
      thresholds = signal_superior_thresholds, figure="HYPER", sample_detail = local_sample_detail,
       probe_features = probe_features)

    hypo_result <- semseeker:::analyze_single_sample( values = signal_values,
      thresholds = signal_inferior_thresholds, figure="HYPO", sample_detail = local_sample_detail,
      probe_features = probe_features)

    both_result_mutations <- semseeker:::analyze_single_sample_both( sample_detail =  local_sample_detail, "MUTATIONS")

    both_result_lesions <- semseeker:::analyze_single_sample_both( sample_detail =  local_sample_detail, "LESIONS")

    delta_result <- delta_single_sample ( values = signal_values, high_thresholds = signal_superior_thresholds,
      low_thresholds = signal_inferior_thresholds, sample_detail = local_sample_detail,
      signal_medians = signal_median_values, probe_features = probe_features)

    deltar_result <- deltar_single_sample ( values = signal_values, high_thresholds = signal_superior_thresholds,
      low_thresholds = signal_inferior_thresholds, sample_detail = local_sample_detail,
      probe_features = probe_features)

    # summary signal
    names(signal_values) <- row.names(signal_data)
    signal_sample <- signal_single_sample( signal_values,local_sample_detail,probe_features)


    sample_status_temp <- dplyr::bind_cols( hyper_result, hypo_result,
      both_result_mutations,both_result_lesions, deltar_result, signal_sample)
    sample_status_temp$Sample_ID <- local_sample_detail$Sample_ID
    # count the real values available for the sample
    sample_status_temp$PROBES_COUNT <- length(stats::na.omit(signal_values))

    if(ssEnv$showprogress)
      progress_bar(sprintf("sample: %s",local_sample_detail$Sample_ID))
    sample_status_temp
  }

  # progress_bar$terminate()
  summary_population <- as.matrix.data.frame(summary_population)
  summary_population <- as.data.frame(summary_population)
  rownames(summary_population) <- summary_population$Sample_ID
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Row count result:", nrow(summary_population))
  rm(signal_data)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population analysis ")
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population with summary - Time taken: ", round(time_taken), " minutes.")

  return(summary_population)
}


