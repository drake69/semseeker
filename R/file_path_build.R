file_path_build <- function(baseFolder, detailsFilename, extension, add_gz = FALSE){

  detailsFilename <- as.vector(sapply(detailsFilename, as.character))

  detailsFilename <- paste0(detailsFilename, collapse="_")

  detailsFilename <- name_cleaning(detailsFilename)

  if(extension!="")
    fileName <- paste0( detailsFilename,".",extension, sep="")

  if(add_gz){
    fileName <- paste0(fileName, ".gz")
  }

  # replace double dots with single dot
  fileName <- gsub("\\.\\.", ".", fileName)

 fp <-  file.path(baseFolder, fileName)

 return(fp)

}



