update_session_info <- function(ssEnv)
{
  # .pkgglobalenv
  assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  #save principal copy
  if (!dir.exists(file.path(ssEnv$session_folder)))
    stop("I'm STOPPING HERE! ", ssEnv$result_folderData)
  else
    saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))
}
