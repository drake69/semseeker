log_event <- function(...)
{
  ssEnv <- get_session_info()
  # append log_event to log file
  log_events <- list(...)
  log_event_to_save <- ""
  for (i in 1:length(log_events))
  {
    log_event_to_save <- paste0(log_event_to_save, log_events[i], sep =" ")
  }

  log_event_to_save <- gsub("  "," ", log_event_to_save)
  file_name <- paste(as.character(Sys.info()["nodename"]),"_session_output.log", sep="")
  log_file <- file.path(ssEnv$session_folder,file_name)

  log_event_to_print <- log_event_to_save
  if (grepl("^BANNER", log_event_to_save))
  {
    log_event_to_print_dash <- "##############################################################################"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",
      gsub("BANNER: ", "", log_event_to_print), "\n",
      log_event_to_print_dash, '\n')
  }

  if (grepl("ERROR:", log_event_to_save))
  {
    log_event_to_print_dash <- "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",  log_event_to_print  , "\n",log_event_to_print_dash)
  }

  if(is.null(ssEnv$verbosity))
    verbosity <- 4
  else
    verbosity <- as.numeric(ssEnv$verbosity)

  # if verbosity is 1, only print ERROR messages
  if (verbosity == 1 && grepl("^ERROR", log_event_to_save) || grepl("^BANNER", log_event_to_save) && !testthat::is_testing())
    message(log_event_to_print)

  # if verbosity is 2, print ERROR and WARNING messages
  if (verbosity == 2 && (grepl("^ERROR", log_event_to_save) || grepl("^WARNING", log_event_to_save)|| grepl("^BANNER", log_event_to_save)) && !testthat::is_testing())
    message(log_event_to_print)

  # if verbosity is 3, print ERROR, WARNING and INFO messages
  if (verbosity == 3 && (grepl("^ERROR", log_event_to_save) || grepl("^WARNING", log_event_to_save) || grepl("^INFO", log_event_to_save) || grepl("^BANNER", log_event_to_save)) && !testthat::is_testing())
    message(log_event_to_print)

  # if verbosity is 4, print ERROR, WARNING, INFO and DEBUG messages
  if (verbosity == 4 && !testthat::is_testing())
    message(log_event_to_print)

  cat(log_event_to_save, "\n", file = log_file, append = TRUE)
}
