close_env <- function()
{

  ssEnv <- get_session_info()

  if (ssEnv$showprogress)
    progressr::handlers()
  # progressr::handlers(global = FALSE)
  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  message("INFO: ", Sys.time(), " Job Completed !")

}
