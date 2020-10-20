

#' Title
#'
#' @param dataToDump
#' @param fileExtension
#' @param resultFolder
#' @param resultSubFolder
#' @param sampleName
#' @param multipleFileColNames
#'
#' @return
#' @export
#'
#' @examples
dumpSampleAsBedFile <- function(dataToDump, fileExtension, resultFolder, resultSubFolder, sampleName, multipleFileColNames) {


  if (resultSubFolder != "" && !dir.exists(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))) {
    dir.create(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))
  }

  if (!empty(dataToDump) && !startsWith(x = as.character(dataToDump[1, "CHR"]), prefix = "CHR")) {
    chr <- rep(x = "chr", dim(dataToDump)[1])
    chr <- paste0(chr, dataToDump[, "CHR"], sep = "")
    dataToDump[, "CHR"] <- chr
  }

  filePath <- paste0(resultFolder, "/", resultSubFolder, "/", sampleName, fileExtension, sep = "")
  write.table(dataToDump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  if (!empty(dataToDump)) {
    message(sampleName, " ", "saved file ", filePath, " ", Sys.time())

    oldFilePath <- paste0(resultFolder, "/", resultSubFolder, "/", "MULTIPLE", fileExtension, sep = "")

    ## needed random extension to avoid concurrency problem to access same file during parallel operations
    tempExtension <- stri_rand_strings(1, 10)
    filePath <- paste0(resultFolder, "/", resultSubFolder, "/", "MULTIPLE", fileExtension, tempExtension, sep = "")
    sampleNames <- rep(sampleName, dim(dataToDump)[1])
    dataToDump <- data.frame(dataToDump, sampleNames)

    colnames(dataToDump) <- multipleFileColNames

    resultMultiple <- dataToDump
    if (!file.exists(filePath)) {
      resultMultiple <- dataToDump
    } else {
      # oldFile <- read.table(file = filePath,sep = '\t', header = FALSE, quote = '') colnames(oldFile) <- multipleFileColNames resultMultiple <- rbind(oldFile, dataToDump)
    }

    write.table(resultMultiple, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

    ## append to multiple file the single sample file though system, afterwards remove it!
    system(paste0("cat ", filePath, " >> ", oldFilePath, sep = ""))
    system(paste0("rm ", filePath, sep = ""))
    message(sampleName, " ", "saved files ", fileExtension, Sys.time())
  }

}

