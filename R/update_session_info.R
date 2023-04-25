update_session_info <- function(ssEnv)
{
  assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  ssEnv <- .pkgglobalenv$ssEnv
  if(length(ssEnv)>0 & length(ssEnv$session_folder)>0)
    saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))

}
