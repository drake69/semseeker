signal_save <- function(signal_data, sample_sheet, batch_id )
{
  ssEnv <- get_session_info()

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data,")

  # save signal as rds and as pivot
  pivot_file_name <- pivot_file_name("SIGNAL","MEAN","PROBE","")
  readr::write_delim(signal_data, pivot_file_name, col_names = T, progress = FALSE)

  saveRDS(signal_data,file_path_build(ssEnv$result_folderData, c(batch_id,"_signal_data"),"rds"))
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data,")

}
