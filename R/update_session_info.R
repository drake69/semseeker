update_session_info <- function(ssEnv)
{
  assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  if(length(ssEnv)>0 & length(ssEnv$session_folder)>0)
  {
    #save principal copy
    saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))
    # result_folder <- getwd()

    #save backup copy
    # session_folder <- file.path(result_folder,"Log")
    # saveRDS(ssEnv, file.path(result_folder,"session_info.rds"))
  }
}
