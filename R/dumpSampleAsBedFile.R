#' given data and colnames dump as bed file
#'
#' @param dataToDump data frame to dump into bed file with CHR, START, END
#' @param fileName name of the file to save data in
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
  dataToDump$START <- as.numeric(dataToDump$START) - 1
  dataToDump$END <- as.numeric(dataToDump$END)
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


