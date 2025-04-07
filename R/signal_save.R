signal_save <- function(signal_data, sample_sheet, batch_id )
{
  ssEnv <- get_session_info()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data.")
  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION","WHOLE")
  if(file.exists(pivot_file_name))
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data already saved.")
    return()
  }

  signal_data <- signal_data[,unique(sample_sheet$Sample_ID)]

  signal_data$AREA <- rownames(signal_data)

  # move AREA to the first column
  signal_data <- signal_data[, c(ncol(signal_data), 1:(ncol(signal_data)-1))]

  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE","WHOLE")
  polars::as_polars_df(signal_data)$write_parquet(pivot_file_name)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data saved with probe.")

  rm(signal_data)

  signal_data <- polars::pl$scan_parquet(pivot_file_name)
  pp <- polars::as_polars_df(probe_features_get("PROBE"))$lazy()
  # rename AREA as PROBE
  signal_data <- signal_data$with_columns(polars::pl$col("AREA")$alias("PROBE"))
  signal_data <- pp$join(
    signal_data,
    on = c("PROBE"),
    how = "inner"
  )
  signal_data <- signal_data$drop(c("PROBE","PROBE_WHOLE","AREA"))
  signal_data <- signal_data$sort(c("CHR", "START","END"), descending = c(FALSE, FALSE,FALSE))
  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION","WHOLE")

  signal_data$collect()$write_parquet(pivot_file_name)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")

  gc()
}
