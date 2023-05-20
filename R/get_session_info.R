get_session_info <- function()
{
  ssEnv <- .pkgglobalenv$ssEnv
  if(is.null(.pkgglobalenv) | is.null(ssEnv)){
    result_folder <- getwd()
    session_folder <- file.path(result_folder,"Log")
    ssEnv <- list()
    if(file.exists(file.path(session_folder,"session_info.rds")))
      ssEnv <- readRDS( file.path(session_folder,"session_info.rds"))
    assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  }
  return(ssEnv)
}
