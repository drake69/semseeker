get_session_info <- function(result_folder=NULL)
{
  ssEnv <- .pkgglobalenv$ssEnv
  if(is.null(.pkgglobalenv) | is.null(ssEnv)){
    if(is.null(result_folder))
      stop("ERROR: get_session_info called without result folder!")
    session_folder <- dir_check_and_create(result_folder,"Log")
    ssEnv <- list()
    if(file.exists(file.path(session_folder,"session_info.rds")))
      ssEnv <- readRDS( file.path(session_folder,"session_info.rds"))
    ssEnv$session_folder <- session_folder
    assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  }
  return(ssEnv)
}
