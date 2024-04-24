file_path_build <- function(baseFolder, detailsFilename, extension, add_gz = FALSE){

  detailsFilename <- as.vector(sapply(detailsFilename, as.character))
  # remove empty strings
  detailsFilename <- detailsFilename[detailsFilename != ""]
  detailsFilename <- detailsFilename[detailsFilename != " "]
  detailsFilename <- paste0(detailsFilename, collapse = "_")
  # remove final _
  detailsFilename <- gsub("_$","", detailsFilename)
  fileName <- paste0( detailsFilename,".",extension, sep="")
  fileName <- gsub("__","_", fileName)

  if(add_gz){
    fileName <- paste0(fileName, ".gz")
  }

  # remove double dots
  fileName <- gsub("\\.\\.",".", fileName)

  file.path(baseFolder, fileName)

}



