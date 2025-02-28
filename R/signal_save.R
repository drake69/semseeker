signal_save <- function(signal_data, sample_sheet, batch_id )
{
  ssEnv <- get_session_info()

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data.")

  signal_data <- signal_data[,unique(sample_sheet$Sample_ID)]

  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE","")
  polars::as_polars_df(signal_data)$write_parquet(pivot_file_name)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")

  gc()
}
