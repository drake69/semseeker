#' Title
#'
#' @param fileExtension
#' @param resultFolder
#' @param multipleFileColNames
#'
#' @return
#' @export
#'
#' @examples
mergeMultipleBed <- function(populations,figures,anomalies, fileExtension, resultFolder,  multipleFileColNames) {

  for (pop in populations) {
    for (fig in figures) {
      for (anomal in anomalies) {

        souceFolder <- paste(resultFolder,"/", pop, "/", anomal, "_", fig, "/", sep = "")

        sourceFiles <- paste0(souceFolder, "/*",  fileExtension,".temp", sep = "")
        destinationFile <- paste0(souceFolder, "/","MULTIPLE", ".", fig, ".", anomal,".bed", sep = "")

        system(paste0("echo '" ,paste(multipleFileColNames, collapse = "\t") ,  "'  > ", destinationFile, sep = ""))
        system(paste0("cat ", sourceFiles, " >> ", destinationFile, sep = ""))
        system(paste0("rm ", sourceFiles, sep = ""))
      }
    }
  }

}
