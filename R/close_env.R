close_env <- function()
{

  progressr::handlers()
  progressr::handlers(global = FALSE)
  future::plan( future::sequential)
  message("INFO: ", Sys.time(), " Job Completed !")
}
