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

  ssEnv <- get_session_info()
  # #
  start_time <- Sys.time()
  log_event("INFO: ", Sys.time(), " AnalyzePopulation warmingUP ")

  signal_data <- stats::na.omit(signal_data)

  ### get signal_values ########################################################
  sample_sheet <- sample_sheet[order(sample_sheet[, "Sample_ID"], decreasing = FALSE), ]
  existent_samples <- colnames(signal_data)
  sample_names <- sample_sheet$Sample_ID
  missed_samples <- setdiff(setdiff(sample_names, existent_samples), "PROBE")

  if (length(missed_samples) != 0) {
    log_event("INFO: ", Sys.time(), " These samples data are missed: ", paste0(missed_samples, sep = " "))
  }

  log_event("INFO: ", Sys.time(), " WarmedUP AnalyzePopulation")
  log_event("INFO: ", Sys.time(), " Start population analysis")

  # progress_bar <- progress::progress_bar$new(
  #   format = paste("INFO: Performing population analysis [:bar] :percent eta: :eta"),
  #   total = nrow(sample_sheet),
  #   clear = FALSE,
  #   width= 60)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  variables_to_export <- c("sample_sheet", "signal_data", "analyze_single_sample", "ssEnv",
                            "signal_superior_thresholds","deltar_single_sample","signal_inferior_thresholds","iqr","signal_median_values",
                           "bt","bonferroni_threshold", "probe_features", "analyze_single_sample_both", "delta_single_sample", "progress_bar",
                           "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_single_sample")
  i <- 1

  if (!ssEnv$signal_intrasample)
  {
    signal_superior_thresholds <- signal_thresholds$signal_superior_thresholds
    signal_inferior_thresholds <- signal_thresholds$signal_inferior_thresholds
    iqr <- signal_thresholds$iqr
    signal_median_values <- signal_thresholds$signal_median_values
  }

  # for(i in 1:nrow(sample_sheet)) {
  summary_population <-  foreach::foreach(i =1:nrow(sample_sheet), .combine= "rbind", .export = variables_to_export) %dorng% {
    local_sample_detail <- sample_sheet[i,]

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

    hyper_result <- analyze_single_sample( values = signal_values,
      thresholds = signal_superior_thresholds, figure="HYPER", sample_detail = local_sample_detail,
       probe_features = probe_features)

    hypo_result <- analyze_single_sample( values = signal_values,
      thresholds = signal_inferior_thresholds, figure="HYPO", sample_detail = local_sample_detail,
      probe_features = probe_features)

    both_result_mutations <- analyze_single_sample_both( sample_detail =  local_sample_detail, "MUTATIONS")

    both_result_lesions <- analyze_single_sample_both( sample_detail =  local_sample_detail, "LESIONS")

    delta_result <- delta_single_sample ( values = signal_values, high_thresholds = signal_superior_thresholds,
      low_thresholds = signal_inferior_thresholds, sample_detail = local_sample_detail,
      signal_medians = signal_median_values, probe_features = probe_features)

    deltar_result <- deltar_single_sample ( values = signal_values, high_thresholds = signal_superior_thresholds,
      low_thresholds = signal_inferior_thresholds, sample_detail = local_sample_detail,
      probe_features = probe_features)

    # summary signal
    names(signal_values) <- row.names(signal_data)
    signal_sample <- signal_single_sample( signal_values,local_sample_detail,probe_features)

    sample_status_temp <- c( "Sample_ID"=local_sample_detail$Sample_ID, delta_result, hyper_result, hypo_result,
      "MUTATIONS_BOTH"=both_result_mutations,"LESIONS_BOTH"=both_result_lesions, deltar_result, signal_sample)

    if(ssEnv$showprogress)
      progress_bar(sprintf("sample: %s",local_sample_detail$Sample_ID))
    sample_status_temp
  }

  # progress_bar$terminate()
  summary_population <- as.matrix.data.frame(summary_population)
  summary_population <- as.data.frame(summary_population)
  colnames(summary_population) <- c("Sample_ID","DELTAS_HYPO","DELTAS_HYPER","DELTAS_BOTH","MUTATIONS_HYPER","LESIONS_HYPER","PROBES_COUNT","MUTATIONS_HYPO",
                                    "LESIONS_HYPO","MUTATIONS_BOTH","LESIONS_BOTH","DELTAR_HYPO","DELTAR_HYPER","DELTAR_BOTH","BETA_MEAN")
  rownames(summary_population) <- summary_population$Sample_ID
  log_event("INFO: ", Sys.time(), " Row count result:", nrow(summary_population))
  rm(signal_data)

  log_event("INFO: ", Sys.time(), " Completed population analysis ")
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  log_event("INFO: ", Sys.time(), " Completed population with summary - Time taken: ", time_taken)

  return(summary_population)
}


