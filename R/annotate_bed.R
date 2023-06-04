#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param populations vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param anomalies vector of lesions/mutations to use to build the folder path
#' @param groups vector of genomic area to cycle and group the annotated data
#' @param probes_prefix prefix to use to get the annotated probes dataset
#' @param columnLabel label of the column of the genomic area gene, island ,dmr etc..
#' @param groupingColumnLabel label of the column of the genomic sub area body, tss1500
#'
#' @return original bed with genomic area infos
#' @importFrom doRNG %dorng%

annotate_bed <- function (
    populations,
    figures,
    anomalies,
    groups,
    probes_prefix,
    columnLabel,
    groupingColumnLabel)
{

  ssEnv <- get_session_info()

  i <- 0
  final_bed <- NULL
  bedFileName <- file_path_build(ssEnv$result_folderData , c(columnLabel, "ANNOTATED"),"fst")

  if("BETA" %in% anomalies)
    figures <- c(figures,"MEAN")

  anomalies <- anomalies[(anomalies !="BETA")]
  figures <- figures[(figures!="MEAN")]

  # if(file.exists(bedFileName))
  # {
  #   if(file.info(bedFileName)$size < 10)
  #   {
  #     final_bed <- NULL
  #     message("WARNING: ", Sys.time(), " Given up file:", final_bed, " is empty!")
  #   }
  #   else
  #   {
  #     # final_bed <- utils::read.table(bedFileName, stringsAsFactors = TRUE, sep="\t", header = TRUE)
  #     final_bed <- fst::fst(path = bedFileName)
  #     final_bed$VALUE = as.numeric(final_bed$VALUE)
  #   }
  #
  #   final_bed_temp <- unique(final_bed[final_bed$ANOMALY %in% anomalies & final_bed$FIGURE %in% figures,])
  #   if(!plyr::empty(final_bed_temp)
  #     & any(!(anomalies %in% unique(final_bed$ANOMALY)))
  #     & any(!(figures %in% unique(final_bed$FIGURE)))
  #     )
  #     return(final_bed_temp)
  #   else
  #     final_bed_temp <- as.data.frame(final_bed)
  # }

  ssEnv$keysLocal <-
    expand.grid(
      "POPULATION" = populations,
      "FIGURE" = figures,
      "ANOMALY" = anomalies,
      "GROUP" = groups
    )

  message("INFO: ", Sys.time(), " Annotating genomic areas.")
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(ssEnv$keysLocal))

  variables_to_export <- c("ssEnv", "probes_prefix", "dir_check_and_create", "read_multiple_bed", "columnLabel",
                           "groupingColumnLabel", "progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","probes_get")

  # for(i in 1:nrow(ssEnv$keysLocal))
  final_bed <- foreach::foreach(i=1:nrow(ssEnv$keysLocal), .combine = rbind, .export = variables_to_export) %dorng%
    {
      anomal <- as.character(ssEnv$keysLocal[i,"ANOMALY"])
      pop <- as.character(ssEnv$keysLocal[i,"POPULATION"])
      fig <- as.character(ssEnv$keysLocal[i,"FIGURE"])
      grp <- as.character(ssEnv$keysLocal[i,"GROUP"])

      probes <- probes_get(probes_prefix, grp)

      resFolder <- dir_check_and_create(ssEnv$result_folderData,pop)
      tempFile <- read_multiple_bed(anomalyLabel =  anomal, figureLable =  fig, probe_features =  probes,
                                    columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
      if(ssEnv$showprogress)
        progress_bar()
      tempFile
    }

  # colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","PROBE"))
  colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","k27","k450","k850"))
  final_bed <- unique(final_bed[, colname_to_preserve])

  # if(exists("final_bed_temp"))
  #   if(!plyr::empty(final_bed_temp))
  #   {
  #     final_bed <- unique(rbind(final_bed, final_bed_temp))
  #   }

  if(exists("final_bed"))
  {
    final_bed <- as.data.frame(final_bed)
    if(!plyr::empty(final_bed))
      fst::write_fst( x = final_bed,path =  bedFileName, compress = T )
  }

  # utils::write.table(final_bed,bedFileName, row.names = FALSE, sep = "\t", col.names = TRUE)

  # final_bed <- unique(final_bed[, c("SAMPLEID","PROBE",columnLabel,groupingColumnLabel,"FIGURE","ANOMALY","POPULATION","VALUE")])

  return(final_bed)

}

