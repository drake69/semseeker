#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param envir semseekere working infos
#' @param populations vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param anomalies vector of lesions/mutations to use to build the folder path
#' @param groups vector of genomic are to cycle and group the annotated data
#' @param probes_prefix prefix to use to get the annotated probes dataset
#' @param columnLabel label of the column of the genomic area gene, island ,dmr etc..
#' @param groupingColumnLabel label of the column of the genomic sub area body, tss1500
#'
#' @return original bed with genomic area infos
#' @importFrom foreach %dopar%

annotate_bed <- function (
  envir,
  populations ,
  figures ,
  anomalies ,
  groups ,
  probes_prefix ,
  columnLabel ,
  groupingColumnLabel)
  {

  i <- 0
  final_bed <- NULL
  bedFileName <- file_path_build(envir$result_folderData , c(columnLabel, "ANNOTATED"),"csv")

  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
      {
        final_bed <- NULL
        message("Given up file:", final_bed, " is empty!")
      }
    else
      {
        final_bed <-    utils::read.table(bedFileName, stringsAsFactors = TRUE, sep="\t", header = TRUE)
        final_bed$VALUE = as.numeric(final_bed$VALUE)
      }

    final_bed <- final_bed[final_bed$ANOMALY %in% anomalies & final_bed$FIGURE %in% figures,]
    if(!plyr::empty(final_bed))
      return(final_bed)
  }

  envir$keysLocal <-
    expand.grid(
      "POPULATION" = populations,
      "FIGURE" = figures,
      "ANOMALY" = anomalies,
      "GROUP" = groups
    )

  variables_to_export <- c("envir", "probes_prefix", "dir_check_and_create", "read_multiple_bed", "columnLabel", "groupingColumnLabel")

  # for(i in 1:nrow(envir$keysLocal))
  final_bed <- foreach::foreach(i=1:nrow(envir$keysLocal), .combine = rbind, .export = variables_to_export) %dopar%
  {
    anomal <- envir$keysLocal[i,"ANOMALY"]
    pop <- envir$keysLocal[i,"POPULATION"]
    fig <- envir$keysLocal[i,"FIGURE"]
    grp <- envir$keysLocal[i,"GROUP"]

    probes <- get(paste0(probes_prefix, grp,sep=""))
    resFolder <- dir_check_and_create(envir$result_folderData,pop)
    tempFile <- read_multiple_bed(envir=envir, anomalyLabel =  anomal, figureLable =  fig, probe_features =  probes,
                                  columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
    tempFile
  }

  colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","PROBE"))
  final_bed <- final_bed[, colname_to_preserve]
  utils::write.table(final_bed,bedFileName, row.names = FALSE, sep = "\t", col.names = TRUE)

  gc()
  return(final_bed)

}

