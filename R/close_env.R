close_env <- function(envir)
{

  if (envir$showprogress)
    progressr::handlers()
  # progressr::handlers(global = FALSE)
  future::plan( future::sequential)
  message("INFO: ", Sys.time(), " Job Completed !")
}
