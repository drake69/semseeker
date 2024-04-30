update_session_info <- function(ssEnv)
{
  # .pkgglobalenv

  # get ssEnv in memory and overwrite only if lenght is ugual or greater
  ssEnv_local <- semseeker:::get_session_info()
  if (length(ssEnv)>=length(ssEnv_local))
    {
    assign("ssEnv", ssEnv, envir=.pkgglobalenv)
    # options("ssEnv"=ssEnv)
  }


  #save principal copy
  if (!dir.exists(file.path(ssEnv$session_folder)))
    stop("I'm STOPPING HERE! ", ssEnv$result_folderData)
  else
  {
    saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))
    # format date as YYYYMMDD
    today <- Sys.Date()
    saveRDS(ssEnv, file.path(ssEnv$session_folder,paste0(today, "_session_info.rds", sep="")))
  }
}
