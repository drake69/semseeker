#' createPivotResultFromMultipleBed load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param anomalyLabel anomaly definition used to lable folder and files eg
#' MUTATIONS, LESIONS
#' @param figureLable figure used to create the file HYPO HYPER
#' @param probeFeatures features of probe CHR and START and NAME
#'
#' @return list of pivot by CHR and by Gene

#'
createPivotResultFromMultipleBed <- function(anomalyLabel, figureLable, probeFeatures) {

  # TODO: check sample name is a column of the bedfile

  souceFolder <- paste0(ssEnv$resultFolderData, "/", anomalyLabel, "_", figureLable, "/", sep = "")
  fileName <- paste0(souceFolder, "/", "MULTIPLE", ".", figureLable, ".", anomalyLabel, ".bed", sep = "")
  sourceData <- utils::read.table(fileName, sep = "\t", blank.lines.skip = TRUE, fill = TRUE)
  colnames(sourceData) <- c("CHR", "START", "END", "SAMPLEID")

  sourceData$CHR <- as.factor(sourceData$CHR)
  probeFeatures$CHR <- as.factor(paste0("chr", probeFeatures$CHR))
  probeFeatures <- probeFeatures[(probeFeatures$CHR %in% unique((sourceData$CHR))), ]
  droplevels(probeFeatures$CHR)
  droplevels(sourceData$CHR)

  sourceDatabyCHR <- plyr::count(df = sourceData, vars = c("SAMPLEID", "CHR"))
  finalResult <- reshape2::dcast(data = sourceDatabyCHR, SAMPLEID ~ CHR, value.var = "freq")
  finalResult[is.na(finalResult)] <- 0
  finalResult[finalResult == ""] <- 0

  resultByCHR <- finalResult
  rm(finalResult)

  sourceData$CHR <- as.factor(sourceData$CHR)
  sourceData$START <- as.integer(sourceData$START)

  sourceData <- dplyr::left_join(sourceData, probeFeatures, by = c("CHR", "START"))
  sourceData[is.na(sourceData)] <- ""
  sourceData$CHR <- paste0(sourceData$CHR, ifelse(test = is.na(sourceData$GENE) | is.null(sourceData$GENE) | sourceData$GENE == "", yes = paste0(":", sourceData$START), no = "_"), sourceData$GENE, sep = "")

  sourceData$CHR <- as.factor(sourceData$CHR)
  sourceData <- plyr::count(df = sourceData, vars = c("CHR", "SAMPLEID"))
  finalResult <- reshape2::dcast(data = sourceData, CHR ~ SAMPLEID, value.var = "freq", sum)
  finalResult[is.na(finalResult)] <- 0
  finalResult[finalResult == ""] <- 0

  resultByGene <- finalResult

  result <- list(byCHR = resultByCHR, byGene = resultByGene)

  return(result)
}
