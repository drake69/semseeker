#' takes a bed and its location (build with the details of population and genomic area)
#' and annotate with detail about genomic area, the output file is only for the arean and specified subarea
#' @param sample_groups vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param markers vector of lesions/mutations to use to build the folder path
#' @param area vector of genomic area to cycle and group the annotated data
#' @param subarea subarea of the area to annotate..
#'
#' @return original bed with genomic area infos
#' @importFrom doRNG %dorng%

annotate_bed <- function (
    sample_groups,
    figures,
    markers,
    area,
    subarea)
{

  ssEnv <- get_session_info()

  # area and subarea are defined using the filename
  i <- 0
  bedFileName <- file_path_build(ssEnv$result_folderData , c(area,subarea, "ANNOTATED"),"fst")

  sample_groups <- data.frame("SAMPLE_GROUP"=sample_groups)
  localKeys <- reshape::expand.grid.df(ssEnv$keys_areas_subareas_markers_figures, sample_groups)
  localKeys <- subset(localKeys,localKeys$MARKER != "BETA")
  localKeys <- subset(localKeys,localKeys$AREA == area)
  localKeys <- subset(localKeys,localKeys$SUBAREA == subarea)
  localKeys <- subset(localKeys,localKeys$FIGURE %in% figures)
  localKeys <- subset(localKeys,localKeys$MARKER %in%  markers)

  if(file.exists(bedFileName))
  {
    if(file.info(bedFileName)$size < 10)
    {
      message("WARNING: ", Sys.time(), " Given up file:", bedFileName, " is empty!")
    }
    else
    {
      final_bed <- fst::fst(path = bedFileName)
      final_bed$VALUE = as.numeric(final_bed$VALUE)
    }

    final_bed_temp <- as.data.frame(final_bed)
    existing_keys <- unique(final_bed_temp[,c("SAMPLE_GROUP","FIGURE","MARKER")])
    existing_keys$AREA <- area
    existing_keys$SUBAREA <- subarea
    existing_keys$pasted <- apply(existing_keys,1, paste , collapse = "_" )
    localKeys$pasted <-apply(localKeys[,c("SAMPLE_GROUP","FIGURE","MARKER","AREA","SUBAREA")],1, paste , collapse = "_" )
    localKeys <- localKeys[ !( localKeys$pasted %in%  existing_keys$pasted ),]
    if(nrow(localKeys)==0)
      return(final_bed_temp)
    rm(final_bed)
  }


  message("INFO: ", Sys.time(), " Annotating genomic area.")
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))

  variables_to_export <- c("ssEnv", "dir_check_and_create", "read_multiple_bed", "subarea",
                            "progress_bar","progression_index", "progression", "progressor_uuid",
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
      dataToAnnotate <- read_multiple_bed(marker = marker ,sample_group =   sample_group, figure = figure)
      if(is.null(dataToAnnotate))
        next

      if(sum(grepl("END", colnames(probe_features)))>0) #dplyr::inner_join
        dataToAnnotate <- merge(dataToAnnotate, probe_features, by = c("CHR", "START","END"))
      else
        dataToAnnotate <- merge(dataToAnnotate, probe_features, by = c("CHR", "START"))

      dataToAnnotate <-subset(dataToAnnotate, !is.na(eval(parse(text=area_subarea))))

      if(nrow(dataToAnnotate)==0)
        next
      sourceDataColNames <- colnames(dataToAnnotate)
      sourceDataColNames[which(sourceDataColNames==area_subarea)] <- "AREA"
      colnames(dataToAnnotate) <- sourceDataColNames

      dataToAnnotate$CHR <- as.factor(dataToAnnotate$CHR)

      probe_features <- probe_features[(probe_features$CHR %in% unique((dataToAnnotate$CHR))), ]
      droplevels(probe_features$CHR)
      droplevels(dataToAnnotate$CHR)

      if(ssEnv$showprogress)
        progress_bar()

      colname_to_preserve <- !(colnames(dataToAnnotate) %in%  c("START","END","K27","K450","K850"))
      dataToAnnotate <- unique(dataToAnnotate[, colname_to_preserve])

      # dataToAnnotate

      if(exists("final_bed"))
      {
        if(length(colnames(final_bed)) != length(colnames(dataToAnnotate)))
          browser()
        final_bed <- rbind(final_bed, dataToAnnotate)
        }
      else
        final_bed <- dataToAnnotate
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

