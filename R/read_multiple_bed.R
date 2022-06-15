#' read multiple bed with annotated data as per input parameter
#'
#' @param envir semseekere working infos
#' @param anomalyLabel anomaly definition used to label folder and files eg MUTATIONS, LESIONS
#' @param probe_features features of probe CHR and START and NAME
#' @param figureLable figures like hypo/hyper to built the data path
#' @param columnLabel name of column in the annotation dataset to select genomic area (gene, island etc..)
#' @param populationName name of the population used to build the data path
#' @param groupingColumnLabel name of the genomic sub area
#'
#' @return list of pivot by column identified with column Label and by Sample

#'
read_multiple_bed <- function(envir, anomalyLabel, figureLable, probe_features, columnLabel, populationName, groupingColumnLabel)
{
  POSITION <- NULL
  # browser()
  f <- paste0(anomalyLabel,"_", figureLable, sep="")
  souceFolder <- dir_check_and_create(envir$result_folderData, c(as.character(populationName),f))

  if(anomalyLabel=="DELTAS")
  {
    file_extension <- "bedgraph"
    col_names <- c("CHR", "START", "END","DELTAS","SAMPLEID")
  }
  else
  {
    file_extension <- "bed"
    col_names <- c("CHR", "START", "END", "SAMPLEID")
  }


  # browser()
  fileName <-file_path_build(souceFolder,c("MULTIPLE",as.character(anomalyLabel),as.character(figureLable)),file_extension)

  if(!file.exists(fileName))
    return(NULL)

  sourceData <- utils::read.table(fileName, sep = "\t", blank.lines.skip = TRUE, fill = FALSE, col.names = col_names, header = FALSE)
  colnames(sourceData) <- col_names

  sourceData$CHR <- as.factor(sourceData$CHR)

  # probe_features <- subset(probe_features, POSITION == 1)
  probe_features$CHR <- as.factor(paste0("chr", probe_features$CHR))

  probe_features <- probe_features[(probe_features$CHR %in% unique((sourceData$CHR))), ]
  droplevels(probe_features$CHR)
  droplevels(sourceData$CHR)

  sourceData <- dplyr::inner_join(sourceData, probe_features, by = c("CHR", "START"))
  sourceData <-subset(sourceData, !is.na(eval(parse(text=columnLabel))))
  sourceData[is.na(sourceData)] <- ""

  if(anomalyLabel!="DELTAS")
    sourceData <- data.frame(sourceData, "DELTAS"=1)

  sourceData <- sourceData[,c("CHR", "START", "END","DELTAS","SAMPLEID")]

  colnames(sourceData) <- c("CHR", "START", "END","VALUE","SAMPLEID")

  # output with column freq
  # message("multiple nrow data", nrow(sourceData))

  if(nrow(sourceData)==0)
    return(NULL)

  sourceData <- data.frame(sourceData,"FIGURE" = figureLable, "ANOMALY" = anomalyLabel, "POPULATION" = populationName)
  sourceData$POPULATION <- as.factor(sourceData$POPULATION)
  sourceData$ANOMALY <- as.factor(sourceData$ANOMALY)
  sourceData$FIGURE <- as.factor(sourceData$FIGURE)
  # message("read multiple nrow data", nrow(sourceData))
  return(sourceData)

}
