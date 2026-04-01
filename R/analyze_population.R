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
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " AnalyzePopulation warmingUP ")

  nrow_before <- nrow(signal_data)
  signal_data <- stats::na.omit(signal_data)
  nrow_after <- nrow(signal_data)
  if (nrow_before != nrow_after) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Removed ", nrow_before - nrow_after, " rows with NA values")
    stop()
  }

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

  variables_to_export <- c("sample_sheet", "signal_data", "analyze_single_sample", "ssEnv",
    "signal_superior_thresholds","deltar_single_sample","signal_inferior_thresholds","iqr","signal_median_values",
    "bt","bonferroni_threshold", "probe_features", "analyze_single_sample_both", "delta_single_sample", "progress_bar",
    "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_single_sample",
    "get_session_info","bed_file_name","signal_thresholds")
  i <- 1

  for(i in 1:nrow(sample_sheet)) {
  # foreach::foreach(i =1:nrow(sample_sheet), .export = variables_to_export) %dorng% {
    local_sample_detail <- sample_sheet[i,]
    ssEnv <- get_session_info()
    signal_values <- signal_data[,local_sample_detail$Sample_ID]
    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL","MEAN")
    if(!file.exists(bed_filename))
      signal_single_sample( signal_values,local_sample_detail,probe_features)
    if(ssEnv$showprogress)
      progress_bar(sprintf("Saving signal of sample: %s",local_sample_detail$Sample_ID))
  }
  gc()

  progress_bar <- NULL
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  rm(signal_data)
  # for(i in 1:nrow(sample_sheet)) {
  foreach::foreach(i =1:nrow(sample_sheet), .export = variables_to_export) %dorng% {

    local_sample_detail <- sample_sheet[i,]
    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL","MEAN")
    signal_values <- utils::read.delim(bed_filename, header = FALSE, sep = "\t")
    colnames(signal_values) <- c("CHR", "START", "END", "VALUE")

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

    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPER")
    if(!file.exists(bed_filename))
      analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPER", sample_detail = local_sample_detail)

    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPO")
    if(!file.exists(bed_filename))
      analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPO", sample_detail = local_sample_detail)

    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAS","HYPO")
    if(!file.exists(bed_filename))
      delta_single_sample( values = signal_values,thresholds = signal_thresholds , sample_detail = local_sample_detail)

    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAR","HYPO")
    if(!file.exists(bed_filename))
      deltar_single_sample ( values = signal_values, thresholds = signal_thresholds,sample_detail = local_sample_detail)

    if(ssEnv$showprogress)
      progress_bar(sprintf("Performed sample: %s",local_sample_detail$Sample_ID))
  }

  gc()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Row count result:", nrow(sample_sheet))
  if(exists("signal_data"))
    rm(signal_data)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population analysis ")
  end_time <- Sys.time()
  time_taken <- difftime(end_time,start_time, units = "mins")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population with summary - Time taken: ", time_taken, " minutes.")

}


