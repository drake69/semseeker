#' signal_single_sample
#'
#' @param values signal vaòues
#' @param sample_detail detais of sample
#' @param probe_features annotation probe
#'
#' @return signal mean
#'
signal_single_sample <- function(values,sample_detail,probe_features)
{
  ssEnv <- get_session_info()

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData, c(sample_detail$Sample_Group ,paste0("SIGNAL","_", "MEAN", sep = "")))
  # save signal values as bed file per each sample
  # probe_features <- probe_features[probe_features$PROBE== row.names(values),]
  # TO DO: check annotation
  # probe_features <- probe_features[1:length(values),]
  signal_values_annotated <- data.frame(as.data.frame(probe_features), "VALUE" = values, row.names = probe_features$PROBE)[, c("CHR", "START", "END","VALUE")]
  dump_sample_as_bed_file(
    data_to_dump = signal_values_annotated,
    fileName = file_path_build(baseFolder =  folder_to_save, detailsFilename =  c(sample_detail$Sample_ID,"SIGNAL","MEAN"), extension = "bedgraph", add_gz=TRUE)
  )

  result <- ""
  result <- result[-1]
  result["BETA_MEAN"] <- mean(signal_values_annotated$VALUE)


  return(result)

}
