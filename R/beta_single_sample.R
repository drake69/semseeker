#' beta_single_sample
#'
#' @param values beta vaòues
#' @param sample_detail detais of sample
#' @param probe_features annotation probe
#'
#' @return beta mean
#'
beta_single_sample <- function(values,sample_detail,probe_features)
{
  ssEnv <- get_session_info()

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData, c(sample_detail$Sample_Group ,paste0("BETA","_", "MEAN", sep = "")))
  # save beta values as bed file per each sample
  beta_values_annotated <- data.frame(as.data.frame(probe_features), "VALUE" = values, row.names = probe_features$PROBE)[, c("CHR", "START", "END","VALUE")]
  dump_sample_as_bed_file(
    data_to_dump = beta_values_annotated,
    fileName = file_path_build(baseFolder =  folder_to_save, detailsFilename =  c(sample_detail$Sample_ID,"BETA","MEAN"), extension = "bedgraph")
  )

  result <- ""
  result <- result[-1]
  result["BETA_MEAN"] <- mean(beta_values_annotated$VALUE)


  return(result)

}
