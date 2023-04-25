get_session_info <- function(result_folder)
{
  session_folder <- file.path(result_folder,"Log")

  ssEnv <- list()
  if(file.exists(file.path(session_folder,"session_info.rds")))
    ssEnv <- readRDS( file.path(session_folder,"session_info.rds"))

  assign("ssEnv", ssEnv, envir=.pkgglobalenv)
}
