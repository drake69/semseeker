update_session_info <- function(ssEnv)
{

  saveCount <- 0
  if (is.null(ssEnv) | length(ssEnv)==0)
    stop("DEBUG: I'm STOPPING HERE! You called update session info without ssEnv!")

  # save to environment
  assign("ssEnv", ssEnv, envir=.pkgglobalenv)
  saveCount <- saveCount + 1

  # save to file
  saveRDS(ssEnv, file.path(ssEnv$session_folder,"session_info.rds"))
  # format date as YYYYMMDD
  today <- Sys.Date()
  saveRDS(ssEnv, file.path(ssEnv$session_folder,paste0(today, "_session_info.rds", sep="")))
  saveCount <- saveCount + 1

  # save for doFuture
  # options("ssEnv"=ssEnv)
  # saveCount <- saveCount + 1

  if (saveCount != 2)
    stop("I'm STOPPING HERE! ", ssEnv$result_folderData)

  return(ssEnv)
}
