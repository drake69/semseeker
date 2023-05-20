#' given data and colnames dump as bed file
#'
#' @param data_to_dump data frame to dump into bed file with CHR, START, END
#' @param fileName name of the file to save data in
#'
#' @return nothing
#'

dump_sample_as_bed_file <- function(data_to_dump, fileName) {

  ssEnv <- get_session_info()
  # message("dump_sample_as_bed_file ssEnv:", length(ssEnv))
  # message("dump_sample_as_bed_file:", ssEnv$result_folderData)

  if (!plyr::empty(data_to_dump) && !startsWith(x = toupper(as.character(data_to_dump[1, "CHR"])), prefix = "CHR")) {
    chr <- rep(x = "chr", dim(data_to_dump)[1])
    chr <- paste0(chr, data_to_dump[, "CHR"], sep = "")
    data_to_dump[, "CHR"] <- chr
  }

  # bed coordinate must start from zero!
  data_to_dump$START <- as.numeric(data_to_dump$START)
  data_to_dump$END <- as.numeric(data_to_dump$END)
  if (!plyr::empty(data_to_dump)) {

    # message("trying to save: ", fileName)

    # save file bed per sample
    utils::write.table(data_to_dump, file = fileName, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

    # save file bed per sample temporary to reuse for aggregated bed file
    # filePath <- paste(fileName,"",".temp")
    # sample_names <- rep(sampleName, dim(data_to_dump)[1])
    # data_to_dump <- data.frame(data_to_dump, sample_names)
    # colnames(data_to_dump) <- multipleFileColNames
    #
    # utils::write.table(data_to_dump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  # message("dump_sample_as_bed_file: ", fileName)
}


