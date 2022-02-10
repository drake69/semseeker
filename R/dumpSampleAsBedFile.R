#' given data and colnames dump as bed file
#'
#' @param dataToDump data frame to dump into bed file with CHR, START, END
#' @param fileExtension extension to give to written file
#' @param resultSubFolder folder under the root result folder to group files
#' @param sampleName name of the sample to create the file for
#' @param multipleFileColNames names to use as header for the files dumped
#'
#' @return nothing

#'

dumpSampleAsBedFile <- function(dataToDump, fileName) {


  if (!plyr::empty(dataToDump) && !startsWith(x = toupper(as.character(dataToDump[1, "CHR"])), prefix = "CHR")) {
    chr <- rep(x = "chr", dim(dataToDump)[1])
    chr <- paste0(chr, dataToDump[, "CHR"], sep = "")
    dataToDump[, "CHR"] <- chr
  }

  # bed coordinate must start from zero!
  dataToDump$START <- dataToDump$START - 1
  if (!plyr::empty(dataToDump)) {

  # save file bed per sample
    utils::write.table(dataToDump, file = fileName, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


  # save file bed per sample temporary to reuse for aggregated bed file
    # filePath <- paste(fileName,"",".temp")
    # sampleNames <- rep(sampleName, dim(dataToDump)[1])
    # dataToDump <- data.frame(dataToDump, sampleNames)
    # colnames(dataToDump) <- multipleFileColNames
    #
    # utils::write.table(dataToDump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}


