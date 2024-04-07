signal_single_sample <- function( values,sample_detail,probe_features)
{
  ssEnv <- get_session_info()

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData, c(sample_detail$Sample_Group ,paste0("SIGNAL","_", "MEAN", sep = "")))
  # save signal values as bed file per each sample
  signal_values_annotated <- data.frame(as.data.frame(probe_features), "VALUE" = values, row.names = probe_features$PROBE)[, c("CHR", "START", "END","VALUE")]
  dump_sample_as_bed_file(
    data_to_dump = signal_values_annotated,
    fileName = file_path_build(baseFolder =  folder_to_save, detailsFilename =  c(sample_detail$Sample_ID,"SIGNAL","MEAN"), extension = "bedgraph")
  )

  result <- ""
  result <- result[-1]
  result["BETA_MEAN"] <- mean(signal_values_annotated$VALUE)

  signal_intrasample <- FALSE
  if (signal_intrasample)
  {
    q <- stats::quantile(values)
    q1 <- as.numeric(q[2])
    q3 <- as.numeric(q[4])
    y_med <- as.numeric(q[3])
    iqr <- stats::IQR(values)
    iqrmult <- 3
    y_sup <- q3 + iqrmult * iqr
    y_inf <- q1 - iqrmult * iqr
    signal_superior_thresholds <- rep(y_sup,length(values))
    signal_inferior_thresholds <- rep(y_inf,length(values))
  }

  return(result)

}
