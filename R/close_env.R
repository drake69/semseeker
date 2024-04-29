close_env <- function()
{

  ssEnv <- semseeker:::get_session_info()

  if (ssEnv$showprogress)
    progressr::handlers()

  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Job Completed !")
  log_event("DEBUG: --------------------------------------------------------------")

}
