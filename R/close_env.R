close_env <- function()
{

  ssEnv <- get_session_info()

  if (ssEnv$showprogress)
    progressr::handlers()

  # build all the folder tree in result_folder

  remove_empty_folders(ssEnv$result_folder)

  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Job Completed !")
  log_event("DEBUG: --------------------------------------------------------------")

}


remove_empty_folders <- function(path) {
  files <- list.files(path, full.names = TRUE)

  for (file in files) {
    if (file.info(file)$isdir) {
      remove_empty_folders(file)  # Controlla ricorsivamente le sottocartelle

      # Dopo il controllo, se la cartella è vuota, la rimuove
      if (length(list.files(file)) == 0) {
        unlink(file, recursive = TRUE)
        # cat("Rimossa cartella vuota:", file, "\n")
      }
    }
  }
}
