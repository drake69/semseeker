file_path_build <- function(baseFolder, detailsFilename, extension){

  detailsFilename <- as.vector(sapply(detailsFilename, as.character))
  fileName <- paste0( paste0(detailsFilename, collapse = "_", sep=""),".",extension, sep="")
  fileName <- gsub("__","_", fileName)

  file.path(baseFolder, fileName)

}



