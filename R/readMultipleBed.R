#' read multiple bed with annotated data as per input parameter
#'
#' @param resultFolder folder as root for bed files organized per anomaly
#' @param anomalyLabel anomaly definition used to label folder and files eg MUTATIONS, LESIONS
#' @param probeFeatures features of probe CHR and START and NAME
#' @param figureLable figures like hypo/hyper to built the data path
#' @param columnLabel name of column in the annotation dataset to select genomic area (gene, island etc..)
#' @param populationName name of the population used to build the data path
#' @param groupingColumnLabel name of the genomic sub area
#'
#' @return list of pivot by column identified with column Label and by Sample

#'
readMultipleBed <- function(resultFolder, anomalyLabel, figureLable, probeFeatures, columnLabel, populationName, groupingColumnLabel)
{
  POSITION  <- NULL
  souceFolder <- paste(resultFolder, "/", anomalyLabel, "_", figureLable, "/","/",  sep = "")
  fileName <- paste(souceFolder, "/", "MULTIPLE", ".", figureLable, ".", anomalyLabel, ".bed", sep = "")
  if(!file.exists(fileName))
    return(NULL)

  # sourceData <- readr::read_delim(fileName, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  sourceData <- utils::read.table(fileName, sep = "\t", blank.lines.skip = TRUE, fill = FALSE, col.names = c("CHR", "START", "END", "SAMPLENAME"), header = TRUE)
  colnames(sourceData) <- c("CHR", "START", "END", "SAMPLENAME")

  sourceData$CHR <- as.factor(sourceData$CHR)

  probeFeatures <- subset(probeFeatures, POSITION == 1)
  probeFeatures$CHR <- as.factor(paste0("chr", probeFeatures$CHR))

  probeFeatures <- probeFeatures[(probeFeatures$CHR %in% unique((sourceData$CHR))), ]
  droplevels(probeFeatures$CHR)
  droplevels(sourceData$CHR)

  sourceData <- dplyr::inner_join(sourceData, probeFeatures, by = c("CHR", "START"))
  sourceData <-subset(sourceData, !is.na(eval(parse(text=columnLabel))))
  sourceData[is.na(sourceData)] <- ""
  sourceData <- plyr::count(df = sourceData, vars = c(columnLabel, "SAMPLENAME", groupingColumnLabel))
  # output with column freq

  # names(sourceData)[names(sourceData) == "freq"] <- figureLable

  if(dim(sourceData)[1]==0)
    return(NULL)
  sourceData <- data.frame(sourceData,"FIGURE" = figureLable, "ANOMALY" = anomalyLabel, "POPULATION" = populationName)
  sourceData$POPULATION <- as.factor(sourceData$POPULATION)
  sourceData$ANOMALY <- as.factor(sourceData$ANOMALY)
  sourceData$FIGURE <- as.factor(sourceData$FIGURE)
  return(sourceData)

}
