beta_single_sample <- function( values,sample_detail,probe_features)
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

  beta_intrasample <- FALSE
  if (beta_intrasample)
  {
    q <- stats::quantile(values)
    q1 <- as.numeric(q[2])
    q3 <- as.numeric(q[4])
    y_med <- as.numeric(q[3])
    iqr <- stats::IQR(values)
    iqrmult <- 3
    y_sup <- q3 + iqrmult * iqr
    y_inf <- q1 - iqrmult * iqr
    beta_superior_thresholds <- rep(y_sup,length(values))
    beta_inferior_thresholds <- rep(y_inf,length(values))
  }

  return(result)

}
