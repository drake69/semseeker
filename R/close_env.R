close_env <- function()
{

  ssEnv <- get_session_info()

  if (ssEnv$showprogress)
    progressr::handlers()

  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  log_event("INFO: ", Sys.time(), " Job Completed !")
  log_event("DEBUG: --------------------------------------------------------------")

}
