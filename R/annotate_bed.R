#' @title create an annotated file for each marker, figure, area and subarea, each file has all the sample_groups used to calculate epimutation
#'
#' @return nothing
#' @importFrom doRNG %dorng%

annotate_bed <- function ()
{
  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
  i <- 0
  dest_folder <- dir_check_and_create(ssEnv$result_folderData,subFolders = c("Annotated"))
  localKeys <-ssEnv$keys_areas_subareas_markers_figures
  message("INFO: ", Sys.time(), " Annotating genomic area.")

  total_cycles <- nrow(localKeys)* length(ssEnv$keys_sample_groups)
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:total_cycles)
  else
    progress_bar <- ""

  variables_to_export <- c("ssEnv", "dir_check_and_create", "read_multiple_bed", "subarea",
                            "progress_bar","progression_index", "progression", "progressor_uuid",
                            "owner_session_uuid", "trace","probe_features_get","dest_folder", "localKeys",
                            "file_path_build")

  # for(i in 1:nrow(localKeys))
  foreach::foreach(i=1:nrow(localKeys), .export = variables_to_export) %dorng%
    {
      # i <- 1
      for (g in 1:length(ssEnv$keys_sample_groups))
      {
        # g <- 1
        marker <- as.character(localKeys[i,"MARKER"])
        sample_group <- as.character(ssEnv$keys_sample_groups[g,"SAMPLE_GROUP"])
        figure <- as.character(localKeys[i,"FIGURE"])
        subarea <- as.character(localKeys[i,"SUBAREA"])
        area <- as.character(localKeys[i,"AREA"])
        area_subarea <- paste(area,"_", subarea, sep="")

        bedFileName <- file_path_build(dest_folder , c(marker, figure, area,subarea, "Annotated"),"fst")
        if (file.exists(bedFileName))
          next

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

        colname_to_preserve <- !(colnames(dataToAnnotate) %in%  c("START","END","K27","K450","K850"))
        dataToAnnotate <- unique(dataToAnnotate[, colname_to_preserve])

        # dataToAnnotate

        if(exists("final_bed"))
          final_bed <- rbind(final_bed, dataToAnnotate)
        else
          final_bed <- dataToAnnotate

        if(ssEnv$showprogress)
          progress_bar(sprintf("Annotating multiple files."))
      }


      if(exists("final_bed"))
      {
        colname_to_preserve <- !(colnames(final_bed) %in%  c("START","END","K27","K450","K850"))
        final_bed <- unique(final_bed[, colname_to_preserve])
        final_bed <- as.data.frame(final_bed)
        if(!plyr::empty(final_bed))
        {
          fst::write_fst( x = final_bed,path =  bedFileName, compress = T )
          message("INFO: annotated file to ", bedFileName)
        }
      }

    }

}

