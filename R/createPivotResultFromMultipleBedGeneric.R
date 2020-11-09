#' createPivotResultFromMultipleBedGeneric load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param resultFolder folder as root for bedfiles organized per anomaly
#' @param anomalyLabel anomaly definition used to lable folder and files eg
#' MUTATIONS, LESIONS
#' @param figureLable figure used to create the file HYPO HYPER
#' @param probeFeatures features of probe CHR and START and NAME
#' @param columnLabel
#'
#' @return list of pivot by column identified with columnLabel and by Sample
#' @export
#'
createPivotResultFromMultipleBedGeneric <- function(resultFolder, anomalyLabel, figureLable, probeFeatures, columnLabel) {

  POSITION <- NULL
  # TODO: check sample name is a column of the bedfile

  souceFolder <- paste(resultFolder, "/", anomalyLabel, "_", figureLable, "/", sep = "")
  fileName <- paste(souceFolder, "/", "MULTIPLE", ".", figureLable, ".", anomalyLabel, ".bed", sep = "")
  sourceData <- utils::read.table(fileName, sep = "\t", blank.lines.skip = TRUE, fill = TRUE)
  colnames(sourceData) <- c("CHR", "START", "END", "SAMPLENAME")

  sourceData$CHR <- as.factor(sourceData$CHR)

  probeFeatures <- subset(probeFeatures, POSITION == 1)
  probeFeatures$CHR <- as.factor(paste0("chr", probeFeatures$CHR))

  probeFeatures <- probeFeatures[(probeFeatures$CHR %in% unique((sourceData$CHR))), ]
  droplevels(probeFeatures$CHR)
  droplevels(sourceData$CHR)

  sourceData <- dplyr::left_join(sourceData, probeFeatures, by = c("CHR", "START"))
  sourceData <-subset(sourceData, !is.na(eval(parse(text=columnLabel))))
  sourceData[is.na(sourceData)] <- ""

  sourceData[,columnLabel] <- as.factor(sourceData[,columnLabel])
  sourceData <- plyr::count(df = sourceData, vars = c(columnLabel, "SAMPLENAME"))
  finalResult <- reshape2::dcast(data = sourceData, eval(parse(text=columnLabel)) ~ SAMPLENAME, value.var = "freq", sum)
  finalResult[is.na(finalResult)] <- 0
  finalResult[finalResult == ""] <- 0

  return(finalResult)
}

