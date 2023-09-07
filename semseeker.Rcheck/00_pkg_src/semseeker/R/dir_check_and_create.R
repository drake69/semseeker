#' dir_check_and_create
#'
#' @param baseFolder folder to look in
#' @param subFolders sub folders to create, complete tree
#'
#' @return full path
#'
dir_check_and_create <- function  (baseFolder, subFolders)
{
  folderSep <- as.character(.Platform$file.sep)
  parts <- unlist(strsplit(as.character(baseFolder), folderSep))
  # do.call(file.path, as.list(c(parts[1:length(parts) - 1], subFolders)))
  subFolders <- as.vector(sapply(subFolders, as.character))
  subFolders <-c(parts[1:length(parts)], subFolders)

  for( y in 1:length(subFolders))
  {
    subFolder <- subFolders[y]
    browser
    if(!exists("subF"))
      subF <- file.path(subFolder)
    else
    {
      subF <- file.path(subF, subFolder)
      if(!dir.exists(subF))
        dir.create(subF)
    }
  }
  return(file.path(as.character(subF)))
}
