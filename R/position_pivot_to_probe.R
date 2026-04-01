
position_pivot_to_probe <- function(signal_data)
{
  probe_features <- semseeker::pp_tot[,c("PROBE","CHR","START","END","K27","K450","K850")]
  # transform signal_data into parquet
  signal_data_temp <- polars::as_polars_df(signal_data)
  signal_data_temp <- signal_data_temp$lazy()
  cols_to_remove <- colnames(probe_features)
  cols_to_remove <- cols_to_remove[!(cols_to_remove %in% c("PROBE"))]
  probe_features <- polars::as_polars_df(probe_features)
  probe_features <- probe_features$lazy()
  signal_data_temp <- probe_features$join(
    signal_data_temp,
    on = c("CHR", "START", "END"),
    how = "inner"
  )
  signal_data_temp <- signal_data_temp$unique()$collect()
  signal_data_temp <- as.data.frame(signal_data_temp)
  best_tech <- colSums(signal_data_temp[,c("K27","K450","K850")])
  best_tech <-  c("K27","K450","K850")[which(best_tech==max(best_tech))]
  signal_data_temp <- signal_data_temp[signal_data_temp[,best_tech]==TRUE,]
  signal_data_temp <- signal_data_temp[, !colnames(signal_data_temp) %in% cols_to_remove]
  # move PROBE as first column
  signal_data_temp <- signal_data_temp[,c("PROBE",colnames(signal_data_temp)[!(colnames(signal_data_temp)=="PROBE")])]
  signal_data <- as.data.frame(signal_data_temp)
  rownames(signal_data) <- signal_data$PROBE
  # remove PROBE column
  signal_data$PROBE <- NULL
  return(signal_data)
}
