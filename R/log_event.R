log_event <- function(...)
{
  ssEnv <- get_session_info()
  # append log_event to log file
  log_events <- list(...)
  log_event_to_save <- ""
  for (i in 1:length(log_events))
  {
    log_event_to_save <- paste0(log_event_to_save, log_events[i])
  }
  log_file <- file.path(ssEnv$session_folder,"session_output.log")
  if (!grepl("DEBUG", log_event_to_save))
    message(log_event_to_save)

  cat(log_event_to_save, "\n", file = log_file, append = TRUE)
}
