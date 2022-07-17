#' delta_single_sample
#'
#' @param values values of methylation
#' @param high_thresholds highest threshold to use for comparison
#' @param low_thresholds lowest threshold to use for comparison
#' @param sample_detail details of sample to analyze
#' @param beta_medians median to use for calculation
#' @param probe_features genomic position of probes
#' @param envir environment to get globals
#' @return summary detail about the analysis
#'
delta_single_sample <- function (envir, values, high_thresholds, low_thresholds, sample_detail, beta_medians, probe_features) {

  # if (!test_match_order(row.names(values), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and Values!", Sys.time())
  # }
  #
  # if (!test_match_order(row.names(high_thresholds), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and high_thresholds!", Sys.time())
  # }
  #
  # if (!test_match_order(row.names(low_thresholds), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and low_thresholds!", Sys.time())
  # }

  # annotated_data <-   as.data.frame(probe_features)
  #
  # annotated_data$values <- values
  # annotated_data$high_thresholds <- high_thresholds
  # annotated_data$low_thresholds <- low_thresholds
  #
  # annotated_data <- subset(annotated_data, annotated_data$values > annotated_data$high_thresholds)
  # annotated_data$DELTA <- annotated_data$high_thresholds - annotated_data$values

  ### get deltas HYPER #########################################################
  deltas_hyper <- data.frame("DELTA"= values - high_thresholds, row.names = probe_features$PROBE)
  # row.names(deltas_hyper) <- row.names(values)
  colnames(deltas_hyper) <- "DELTA"

  # if (!test_match_order(row.names(deltas_hyper), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes on deltas!", Sys.time())
  # }

  deltasAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltas_hyper)
  # deltasAnnotated_hyperSorted <- deltasAnnotated_hyper
  deltasAnnotated_hyperSorted <- sort_by_chr_and_start(deltasAnnotated_hyper)
  deltasAnnotated_hyperSorted <- subset(deltasAnnotated_hyperSorted, deltasAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPER"),"bedgraph"))

  ### get deltas_hypo HYPER #########################################################
  deltas_hypo <- data.frame("DELTA"=  low_thresholds - values, row.names = probe_features$PROBE)
  colnames(deltas_hypo) <- "DELTA"

  # if (!test_match_order(row.names(deltas_hypo), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes on deltas!", Sys.time())
  # }

  deltasAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltas_hypo, row.names = probe_features$PROBE)
  # deltasAnnotated_hypoSorted <- deltasAnnotated_hypo
  deltasAnnotated_hypoSorted <- sort_by_chr_and_start(deltasAnnotated_hypo)
  deltasAnnotated_hypoSorted <- subset(deltasAnnotated_hypoSorted, deltasAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPO"),"bedgraph"))


  ### get deltas BOTH #########################################################
  deltasAnnotated_both <- rbind(deltasAnnotated_hypoSorted, deltasAnnotated_hyperSorted)
  deltasAnnotated_bothSorted <- sort_by_chr_and_start(deltasAnnotated_both)
  deltasAnnotated_bothSorted <- subset(deltasAnnotated_bothSorted, deltasAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTH"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTH"),"bedgraph"))


  ### get deltas from medians #########################################################

  deltas <- data.frame("DELTA"= round(values - beta_medians,5), row.names = probe_features$PROBE)
  colnames(deltas) <- "DELTA"

  # if (!test_match_order(row.names(deltas), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes on deltas!", Sys.time())
  # }

  # deltasAnnotated <- data.frame(as.data.frame(probe_features), deltas)
  # # deltasAnnotatedSorted <- deltasAnnotated
  # deltasAnnotatedSorted <- sort_by_chr_and_start(deltasAnnotated)
  # deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, deltasAnnotatedSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]

  result <- ""
  result <- result[-1]
  result["DELTAS_HYPO"] <- round(mean(deltasAnnotated_hypoSorted$DELTA),5)
  result["DELTAS_HYPER"] <- round(mean(deltasAnnotated_hyperSorted$DELTA),5)
  result["DELTAS_BOTH"] <- round(mean(deltasAnnotated_bothSorted$DELTA),5)


  return(result)
}






