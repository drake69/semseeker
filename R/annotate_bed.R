#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param sample_groups vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param markers vector of lesions/mutations to use to build the folder path
#' @param area vector of genomic area to cycle and group the annotated data
#' @param probes_prefix prefix to use to get the annotated probe_features dataset
#' @param subarea label of the column of the genomic area gene, island ,dmr etc..
#' @param groupingColumnLabel label of the column of the genomic sub area body, tss1500
#'
#' @return original bed with genomic area infos
#' @importFrom doRNG %dorng%

annotate_bed <- function (
    sample_groups,
    figures,
    markers,
    area)
{

  ssEnv <- get_session_info()

  i <- 0
  final_bed <- NULL
  bedFileName <- file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")

  sample_groups <- data.frame("SAMPLE_GROUP"=sample_groups)
  localKeys <- reshape::expand.grid.df(ssEnv$keys_areas_subareas_markers_figures, sample_groups)
  localKeys <- subset(localKeys,MARKER != "BETA")

  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
    {
      final_bed <- NULL
      message("WARNING: ", Sys.time(), " Given up file:", final_bed, " is empty!")
    }
    else
    {
      final_bed <- fst::fst(path = bedFileName)
      final_bed$VALUE = as.numeric(final_bed$VALUE)
    }

    final_bed_temp <- as.data.frame(final_bed)
    # get existing keys from the existing multiple annotated bed file
    existing_keys <- unique(final_bed_temp[,c("SAMPLE_GROUP","FIGURE","MARKER","AREA","SUBAREA")])
    existing_keys$pasted <- apply(existing_keys,1, paste , collapse = "-" )
    localKeys$pasted <-apply(localKeys,1, paste , collapse = "-" )
    localKeys <- localKeys[ !( localKeys$pasted %in%  existing_keys$pasted ),]
    if(nrow(localKeys)==0)
      return(final_bed_temp)
    rm(final_bed)
  }


  message("INFO: ", Sys.time(), " Annotating genomic area.")
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))

  variables_to_export <- c("ssEnv", "dir_check_and_create", "read_multiple_bed", "subarea",
                            "groupingColumnLabel", "progress_bar","progression_index", "progression", "progressor_uuid",
                            "owner_session_uuid", "trace","probe_features_get")

  for(i in 1:nrow(localKeys))
  # final_bed <- foreach::foreach(i=1:nrow(localKeys), .combine = rbind, .export = variables_to_export) %dorng%
    {
      marker <- as.character(localKeys[i,"MARKER"])
      sample_group <- as.character(localKeys[i,"SAMPLE_GROUP"])
      figure <- as.character(localKeys[i,"FIGURE"])
      subarea <- as.character(localKeys[i,"SUBAREA"])
      area <- as.character(localKeys[i,"AREA"])

      area_subarea <- paste(area,"_", subarea, sep="")
      probe_features <- probe_features_get(area_subarea)
      probe_features$CHR <- as.factor(paste0("chr", probe_features$CHR))

      # annotate file
      resFolder <- dir_check_and_create(ssEnv$result_folderData,sample_group)
      sourceData <- read_multiple_bed(marker = marker ,sample_group =   sample_group, figure = figure, area = area)

      if(sum(grepl("END", colnames(probe_features)))>0) #dplyr::inner_join
        sourceData <- merge(sourceData, probe_features, by = c("CHR", "START","END"))
      else
        sourceData <- merge(sourceData, probe_features, by = c("CHR", "START"))
      sourceData <-subset(sourceData, !is.na(eval(parse(text=area))))

      sourceData$CHR <- as.factor(sourceData$CHR)
      if(nrow(probe_features)==0 | nrow(sourceData)==0)
        return(sourceData)

      probe_features <- probe_features[(probe_features$CHR %in% unique((sourceData$CHR))), ]
      droplevels(probe_features$CHR)
      droplevels(sourceData$CHR)

      if(ssEnv$showprogress)
        progress_bar()

      colname_to_preserve <- !(colnames(tempFile) %in%  c("START","END","K27","K450","K850"))
      tempFile <- unique(tempFile[, colname_to_preserve])

      # tempFile

      if(length(colnames(final_bed)) != length(colnames(tempFile)))
        browser()

      if(exists("final_bed"))
        final_bed <- rbind(final_bed, tempFile)
      else
        final_bed <- tempFile
    }

  # colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","PROBE"))
  colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","K27","K450","K850"))
  final_bed <- unique(final_bed[, colname_to_preserve])

  if(exists("final_bed_temp"))
    if(!plyr::empty(final_bed_temp))
    {
      final_bed <- unique(rbind(final_bed, final_bed_temp))
    }

  if(exists("final_bed"))
  {
    final_bed <- as.data.frame(final_bed)
    if(!plyr::empty(final_bed))
      fst::write_fst( x = final_bed,path =  bedFileName, compress = T )
  }

  # utils::write.table(final_bed,bedFileName, row.names = FALSE, sep = "\t", col.names = TRUE)
  # final_bed <- unique(final_bed[, c("SAMPLEID","PROBE",subarea,groupingColumnLabel,"FIGURE","MARKER","SAMPLE_GROUP","VALUE")])

  return(final_bed)

}

