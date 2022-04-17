#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param envir semseekere working infos
#' @param methylation_data whole matrix of data to analyze.
#' @param sliding_window_size size of the sliding widows to compute epilesions
#' default 11 probes.
#' @param beta_superior_thresholds  data frame to select, from the sample sheet,
#' samples to use as control as study population and as refereces two vectors
#' within the first vector the names of the selection colum and tha second
#' vector the study population selector,
#' @param beta_inferior_thresholds name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param sample_sheet name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param beta_medians name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param envir semseekere working infos
#' @param bonferroni_threshold threshold to define which pValue accept for
#' @param probe_features probes detail from 27 to EPIC illumina dataset
#' lesions definition
#' @return files into the result folder with pivot table and bedgraph.
#' @importFrom foreach %dopar%

analize_population <- function(envir, methylation_data, sliding_window_size, beta_superior_thresholds, beta_inferior_thresholds, sample_sheet, beta_medians, bonferroni_threshold = 0.05, probe_features) {

  #@importFrom foreach %dopar%
  # browser()
  start_time <- Sys.time()
  message("AnalyzePopulation warmingUP ", Sys.time())

  ### get beta_values ########################################################
  sample_sheet <- sample_sheet[order(sample_sheet[, "Sample_ID"], decreasing = FALSE), ]
  existent_samples <- colnames(methylation_data)
  sample_names <- sample_sheet$Sample_ID
  # sampleToSelect <- existent_samples[sample_names %in% existent_samples]
  missed_samples <- setdiff(setdiff(sample_names, existent_samples), "PROBE")

  if (length(missed_samples) != 0) {
    message("These samples data are missed: ", paste0(missed_samples, sep = " "), Sys.time())
  }

  message("WarmedUP AnalyzePopulation", Sys.time())
  message("Start population analyze ", Sys.time())

  # summaryFileName <- file.path(envir$result_folderData, "summary.csv")

  variables_to_export <- c("sample_sheet", "methylation_data", "analyze_single_sample", "envir", "sliding_window_size", "beta_superior_thresholds",
                "bonferroni_threshold", "probe_features", "beta_inferior_thresholds", "analyze_single_sample_both", "delta_single_sample", "beta_medians")
  # , .errorhandling = 'remove'
  # .packages=c("dplyr","semseeker")
  i <- 0
  # summary_population <- foreach::foreach(i = 1:nrow(sample_sheet), .combine="rbind",
  #                                       .export=c(ls(envir=globalenv()),"methylation_data", "sliding_window_size", "beta_superior_thresholds", "bonferroni_threshold", "probe_features", "beta_inferior_thresholds", "beta_medians",
  #                                                 "sample_sheet", "analyze_single_sample", "envir", "analyze_single_sample_both", "delta_single_sample"),
  #                                       .packages=c("dplyr"), .multicombine = FALSE, .errorhandling = "remove") %dopar% {
  summary_population <-  foreach::foreach(i =1:nrow(sample_sheet), .combine='rbind', .export = variables_to_export) %dopar% {
  # for(i in 1:nrow(sample_sheet) ) {
    local_sample_detail <- sample_sheet[i,]
    beta_values <- methylation_data[, local_sample_detail$Sample_ID]
    hyper_result <- analyze_single_sample(envir=envir, values = beta_values, sliding_window_size = sliding_window_size,  thresholds = beta_superior_thresholds, figure="HYPER", sample_detail = local_sample_detail, bonferroni_threshold = bonferroni_threshold, probe_features = probe_features)
    hypo_result <- analyze_single_sample(envir=envir, values = beta_values, sliding_window_size = sliding_window_size,  thresholds = beta_inferior_thresholds, figure="HYPO", sample_detail = local_sample_detail, bonferroni_threshold = bonferroni_threshold, probe_features = probe_features)
    both_result <- analyze_single_sample_both(envir=envir, sample_detail =  local_sample_detail)
    delta_result <- delta_single_sample (envir = envir, values = beta_values, high_thresholds = beta_superior_thresholds, low_thresholds = beta_inferior_thresholds, sample_detail = local_sample_detail, beta_medians = beta_medians, probe_features = probe_features)
    sample_status_temp <- c( "Sample_ID"=local_sample_detail$Sample_ID, delta_result, hyper_result, hypo_result, both_result)
    sample_status_temp <- data.frame(sample_status_temp)
    sample_status_temp <- data.frame(t(sample_status_temp))
    rownames(sample_status_temp) <- c(local_sample_detail$Sample_ID)
    as.data.frame(sample_status_temp)
    # if(exists("summary_population"))
    #   summary_population <- rbind(summary_population, sample_status_temp)
    # else
    #   summary_population <- sample_status_temp
  }

  message("Row count result:", nrow(summary_population))
  rm(methylation_data)
  gc()

  message("Completed population analysis ", Sys.time())
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  message("Completed population with Excel summary", Sys.time(), " Time taken: ", time_taken)

  return(summary_population)
}


