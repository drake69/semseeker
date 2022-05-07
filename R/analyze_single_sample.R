#' analyze_single_sample
#'
#' @param values values of methylation
#' @param sliding_window_size size of window sliding to calculate hypergeometric
#' @param thresholds threshold to use for comparison
#' @param sample_detail details of the sample to analyze
#' @param figure which figure's of sasmple will be analized HYPO or HYPER
#' @param bonferroni_threshold bonferroni threshold to validate pValue
#' @param probe_features probes details to be used
#' @param envir environment to get globals
#' @return list of lesion count  and probes count
#'
#'
analyze_single_sample <- function(envir, values, sliding_window_size, thresholds, figure, sample_detail, bonferroni_threshold = 0.05, probe_features) {

  # browser()
  start_time_single_sample <- Sys.time()
  # message(sample_detail$Sample_ID, " ", "SingleSample Sample analysis warmingUP ", Sys.time())
  result <- ""
  result <- result[-1]

  # colnames(values) <- "VALUE"

  # message(sample_detail$Sample_ID, " ", "SingleSample Sample analysis WarmedUP ...", Sys.time())
  # message(sample_detail$Sample_ID, " ", "SingleSample Start sample analyze ", Sys.time())

  mutation_annotated_sorted <- mutations_get(values, figure,thresholds, probe_features, sample_detail$Sample_ID)
  mutation_annoated_sorted_to_save <- subset(mutation_annotated_sorted, mutation_annotated_sorted$MUTATIONS == 1)[, c("CHR", "START", "END")]

  # message("############# SEARCH")
  # message("############# SEARCH",search())
  # message("############# LS",ls())
  # # browser()
  # message("############# envir$result_folderData:", envir$result_folderData)
  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  dump_sample_as_bed_file(
    data_to_dump = mutation_annoated_sorted_to_save,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS",figure),"bed")
  )
  result[paste("MUTATIONS_", figure, sep="")] <- if (!is.null(mutation_annoated_sorted_to_save)) dim(mutation_annoated_sorted_to_save)[1] else 0

  lesionWeighted <- lesions_get(bonferroni_threshold = bonferroni_threshold, sliding_window_size = sliding_window_size, grouping_column = "CHR", mutation_annotated_sorted = mutation_annotated_sorted)
  folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  dump_sample_as_bed_file(
    data_to_dump = lesionWeighted,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"LESIONS",figure),"bed")
  )

  result[paste("LESIONS_", figure, sep="")] <- if (!is.null(lesionWeighted)) dim(lesionWeighted)[1] else 0
  if(figure=="HYPER")
     result["PROBES_COUNT"] <- dim(probe_features)[1]

  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sample_detail$Sample_ID, " ", "Completed sample ", time_taken)
  return(result)
}



