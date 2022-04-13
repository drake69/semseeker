#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param envir semseekere working infos
#' @param populations vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param anomalies vector of lesions/mutations to use to build the folder path
#' @param groups vector of genomic are to cycle and group the annotated data
#' @param probesPrefix prefix to use to get the annotated probes dataset
#' @param columnLabel label of the column of the genomic area gene, island ,dmr etc..
#' @param groupingColumnLabel label of the column of the genomic sub area body, tss1500
#'
#' @return original bed with genomic area infos
#'

annotateBed <- function (
  envir,
  populations ,
  figures ,
  anomalies ,
  groups ,
  probesPrefix ,
  columnLabel ,
  groupingColumnLabel)
  {

  i <- 0
  finalBed <- NULL
  bedFileName <- file_path_build(envir$resultFolderData , c(columnLabel, "ANNOTATED"),"csv")

  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
      {
        finalBed <- NULL
        message("Given up file:", finalBed, " is empty!")
      }
    else
      {
        finalBed <-    utils::read.table(bedFileName, stringsAsFactors = TRUE, sep="\t", header = TRUE)
        finalBed$freq = as.numeric(finalBed$freq)
      }
    return(finalBed)
  }

  # parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c("readMultipleBed","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body",
  #                                                                                   "PROBES_Gene_TSS1500","PROBES_Gene_TSS200","PROBES_Gene_Whole","PROBES_Gene_ExonBnd","PROBES_Gene_1stExon",
  #                                                                                   "PROBES_DMR_DMR","PROBES_Island_Island","PROBES_Island_N_Shelf","PROBES_Island_S_Shelf","PROBES_Island_N_Shore","PROBES_Island_S_Shore",
  #                                                                                   "PROBES_Island_Whole"))


  envir$keysLocal <-
    expand.grid(
      "POPULATION" = populations,
      "FIGURE" = figures,
      "ANOMALY" = anomalies,
      "GROUP" = groups
    )

  toExport <- c("envir", "probesPrefix", "dir_check_and_create", "readMultipleBed", "columnLabel", "groupingColumnLabel")

  # for(i in 1:nrow(envir$keysLocal))
  finalBed <- foreach::foreach(i=1:nrow(envir$keysLocal), .combine = rbind, .export = toExport) %dopar%
  {
    anomal <- envir$keysLocal[i,"ANOMALY"]
    pop <- envir$keysLocal[i,"POPULATION"]
    fig <- envir$keysLocal[i,"FIGURE"]
    grp <- envir$keysLocal[i,"GROUP"]

    probes <- get(paste0(probesPrefix, grp,sep=""))
    resFolder <- dir_check_and_create(envir$resultFolderData,pop)
    tempFile <- readMultipleBed(envir=envir, anomalyLabel =  anomal, figureLable =  fig, probeFeatures =  probes, columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
    tempFile
    # if(exists("finalBed"))
    #   finalBed <- rbind(finalBed, tempFile)
    # else
    #   finalBed <- tempFile
  }

  # message("Annotated bed:")
  # message(bedFileName)
  # browser()
  utils::write.table(finalBed,bedFileName, row.names = FALSE, sep = "\t", col.names = TRUE)

  # if(nrow(finalBed)==0)
  #  stop("Empty file to annotate !")

  gc()
  return(finalBed)

}

