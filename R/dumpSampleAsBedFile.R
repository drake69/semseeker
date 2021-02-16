#' given data and colnames dump as bed file
#'
#' @param dataToDump data frame to dump into bed file with CHR, START, END
#' @param fileExtension extension to give to written file
#' @param resultFolder root folder to save dumped files
#' @param resultSubFolder folder under the root result folder to group files
#' @param sampleName name of the sample to create the file for
#' @param multipleFileColNames names to use as header for the files dumped
#'
#' @return nothing

#'

dumpSampleAsBedFile <- function(dataToDump, fileExtension, resultFolder, resultSubFolder, sampleName, multipleFileColNames) {

  if (resultSubFolder != "" && !dir.exists(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))) {
    dir.create(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))
  }

  if (!plyr::empty(dataToDump) && !startsWith(x = as.character(dataToDump[1, "CHR"]), prefix = "CHR")) {
    chr <- rep(x = "chr", dim(dataToDump)[1])
    chr <- paste0(chr, dataToDump[, "CHR"], sep = "")
    dataToDump[, "CHR"] <- chr
  }

  if (!plyr::empty(dataToDump)) {

  # save file bed per sample
    filePath <- paste0(resultFolder, "/", resultSubFolder, "/", sampleName, fileExtension, sep = "")
    utils::write.table(dataToDump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

    message(sampleName, " ", "saved file ", filePath, " ", Sys.time())

  # save file bed per sample temporary to reuse for aggregated bed file
    filePath <- paste0(resultFolder, "/", resultSubFolder, "/", sampleName , fileExtension,".temp", sep = "")
    sampleNames <- rep(sampleName, dim(dataToDump)[1])
    dataToDump <- data.frame(dataToDump, sampleNames)
    colnames(dataToDump) <- multipleFileColNames

    utils::write.table(dataToDump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    message(sampleName, " ", "saved files ", fileExtension, Sys.time())
  }
}


