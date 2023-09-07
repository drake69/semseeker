#' analyze_single_sample
#'
#' @param values values of methylation
#' @param sliding_window_size size of window sliding to calculate hypergeometric
#' @param thresholds threshold to use for comparison
#' @param sample_detail details of the sample to analyze
#' @param figure which figure's of sasmple will be analized HYPO or HYPER
#' @param bonferroni_threshold bonferroni threshold to validate pValue
#' @param probe_features probe_features details to be used
#' @return list of lesion count  and probe_features count
#'
#' @importFrom doRNG %dorng%
analyze_single_sample <- function(values, sliding_window_size, thresholds, figure, sample_detail, bonferroni_threshold = 0.05, probe_features) {

  ssEnv <- get_session_info()
  result <- ""
  result <- result[-1]
  # message("analyze_single_sample:", ssEnv$result_folderData)

  mutation_annotated_sorted <- mutations_get(values, figure,thresholds, probe_features, sample_detail$Sample_ID)
  mutation_annotated_sorted_to_save <- subset(mutation_annotated_sorted, mutation_annotated_sorted$MUTATIONS == 1)[, c("CHR", "START", "END")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  # message("analyze_single_sample:",folder_to_save)
  dump_sample_as_bed_file(
    data_to_dump = mutation_annotated_sorted_to_save,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS",figure),"bed")
  )
  result[paste("MUTATIONS_", figure, sep="")] <- if (!is.null(mutation_annotated_sorted_to_save)) dim(mutation_annotated_sorted_to_save)[1] else 0

  lesionWeighted <- lesions_get(bonferroni_threshold = bonferroni_threshold, sliding_window_size = sliding_window_size, grouping_column = "CHR", mutation_annotated_sorted = mutation_annotated_sorted)
  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  dump_sample_as_bed_file(
    data_to_dump = lesionWeighted,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"LESIONS",figure),"bed")
  )

  result[paste("LESIONS_", figure, sep="")] <- if (!is.null(lesionWeighted)) dim(lesionWeighted)[1] else 0
  if(figure=="HYPER")
     result["PROBES_COUNT"] <- dim(probe_features)[1]

  return(result)
}



