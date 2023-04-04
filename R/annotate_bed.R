#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param envir semseekere working infos
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
  bedFileName <- file_path_build(envir$result_folderData , c(columnLabel, "ANNOTATED"),"fst")

  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
    {
      final_bed <- NULL
      message("WARNING: ", Sys.time(), " Given up file:", final_bed, " is empty!")
    }
    else
    {
      # final_bed <- utils::read.table(bedFileName, stringsAsFactors = TRUE, sep="\t", header = TRUE)
      final_bed <- fst::fst(path = bedFileName)
      final_bed$VALUE = as.numeric(final_bed$VALUE)
    }

    final_bed_temp <- final_bed[final_bed$ANOMALY %in% anomalies & final_bed$FIGURE %in% figures,]
    if(!plyr::empty(final_bed_temp))
      return(final_bed_temp)
    else
      final_bed_temp <- as.data.frame(final_bed)
  }

  envir$keysLocal <-
    expand.grid(
      "POPULATION" = populations,
      "FIGURE" = figures,
      "ANOMALY" = anomalies,
      "GROUP" = groups
    )

  message("INFO: ", Sys.time(), "Annotating genomic areas.")
  progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  variables_to_export <- c("envir", "probes_prefix", "dir_check_and_create", "read_multiple_bed", "columnLabel",
                           "groupingColumnLabel", "progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace")

  # for(i in 1:nrow(envir$keysLocal))
  final_bed <- foreach::foreach(i=1:nrow(envir$keysLocal), .combine = rbind, .export = variables_to_export) %dorng%
    {
      anomal <- envir$keysLocal[i,"ANOMALY"]
      pop <- envir$keysLocal[i,"POPULATION"]
      fig <- envir$keysLocal[i,"FIGURE"]
      grp <- envir$keysLocal[i,"GROUP"]

      if(probes_prefix=="PROBES" | probes_prefix=="PROBES_CHR_")
      {
        probes <- semseeker::PROBES_CHR_CHR
        # message("DEBUG: loaded probes: PROBES_CHR_CHR")
      }
      else
      {
        probes_name <- paste0(probes_prefix, grp,sep="")
        probes <- get(probes_name, envir = asNamespace("semseeker"))
        # message("DEBUG: loaded probes:", probes_name)
      }


      resFolder <- dir_check_and_create(envir$result_folderData,pop)
      tempFile <- read_multiple_bed(envir=envir, anomalyLabel =  anomal, figureLable =  fig, probe_features =  probes,
                                    columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
      progress_bar(sprintf("sample=%g",i))
      tempFile
    }

  # colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","PROBE"))
  colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END"))
  final_bed <- final_bed[, colname_to_preserve]

  if(exists("final_bed_temp"))
    if(!plyr::empty(final_bed_temp))
      final_bed <- rbind(final_bed, final_bed_temp)

  if(exists("final_bed"))
  {
    final_bed <- as.data.frame(final_bed)
    if(!plyr::empty(final_bed))
      fst::write_fst( x = final_bed,path =  bedFileName, compress = T )
  }

  # utils::write.table(final_bed,bedFileName, row.names = FALSE, sep = "\t", col.names = TRUE)


  return(final_bed)

}

