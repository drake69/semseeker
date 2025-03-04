signal_save <- function(signal_data, sample_sheet, batch_id )
{
  ssEnv <- get_session_info()

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data.")

  signal_data <- signal_data[,unique(sample_sheet$Sample_ID)]

  signal_data$AREA <- rownames(signal_data)

  # move AREA to the first column
  signal_data <- signal_data[, c(ncol(signal_data), 1:(ncol(signal_data)-1))]

  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE","")
  polars::as_polars_df(signal_data)$write_parquet(pivot_file_name)

  pp <- probe_features_get("PROBE")
  signal_data <- merge(pp,signal_data, by.x="PROBE",by.y="AREA")

  # remove PROBE column
  signal_data <- signal_data[, which(colnames(signal_data) != "PROBE")]
  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION","")
  polars::as_polars_df(signal_data)$write_parquet(pivot_file_name)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")

  gc()
}
