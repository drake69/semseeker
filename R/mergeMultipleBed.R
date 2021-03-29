#' Merge multiple temporary bed into one bed file. Using the input parameters t build the location
#'
#' @param populations vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param anomalies vector of lesions/mutations to use to build the folder path
#' @param resultFolder folder to which build the folder tree and save the annotated bed
#' @param fileExtension extension given to filter the files to merge
#' @param multipleFileColNames
#'
#' @return bed file

#'
#' @examples  populations <- c("Reference","Control","Case")
#' mergeMultipleBed(
#'   populations,
#'   figures = c("METHYLATION"),
#'   anomalies = c("DELTAS"),
#'   fileExtension = ".bedgraph",
#'   resultFolder = resultFolder,
#'   multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME", "VALUE")
#' )

mergeMultipleBed <- function(populations,figures,anomalies, fileExtension, resultFolder,  multipleFileColNames) {

  for (pop in populations) {
    for (fig in figures) {
      for (anomal in anomalies) {

        souceFolder <- paste(resultFolder,"/", pop, "/", anomal, "_", fig, "/", sep = "")

        sourceFiles <- paste0(souceFolder, "/*",  fileExtension,".temp", sep = "")
        destinationFile <- paste0(souceFolder, "/","MULTIPLE", ".", fig, ".", anomal,".bed", sep = "")

        system(paste0("echo '" ,paste(multipleFileColNames, collapse = "\t") ,  "'  > ", destinationFile, sep = ""))

        if (.Platform$OS.type == "windows") {

          command <- paste0("type ", (sourceFiles), " > ",(destinationFile), sep = "")
          command <- gsub ("/","\\\\",command)
          shell(command, intern = TRUE)

          command <- paste0("del ", (sourceFiles), sep = "")
          command <- gsub ("/","\\\\",command)
          #print(command)
          shell(command, intern = TRUE)
          # system2(paste0("type ", shQuote(filePath), " > ",shQuote(summaryFileName), sep = ""))
          # system2(paste0("rm ", filePath, sep = ""))
        } else
        {
          system(paste0("cat ", sourceFiles, " >> ", destinationFile, sep = ""))
          system(paste0("rm ", sourceFiles, sep = ""))
        }
      }
    }
  }

}
