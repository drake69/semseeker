#' @title create an annotated file for each marker, figure, area and subarea, each file has all the sample_groups used to calculate epimutation
#'
#' @return nothing
#' @importFrom doRNG %dorng%

annotate_position_pivots <- function ()
{

  browser()
  start_time <- Sys.time()
  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
  i <- 0
  localKeys <-ssEnv$keys_areas_subareas_markers_figures

  # remove POSITION area
  localKeys <- localKeys[localKeys$AREA != "POSITION",]

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating genomic area.")

  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:(nrow(localKeys)))

  variables_to_export <- c("ssEnv", "dir_check_and_create", "subarea",
    "progress_bar","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","probe_features_get", "localKeys",
    "file_path_build","%>%","get_session_info","log_event")

  # check probe features are avaialable
  probe_features <- probe_features_get("PROBE_WHOLE")


  if(nrow(localKeys)==0)
    return()

  # for(i in 1:nrow(localKeys))
  foreach::foreach(i=1:nrow(localKeys), .export = variables_to_export, .combine = plyr::rbind.fill) %dorng%
    {
      marker <- as.character(localKeys[i,"MARKER"])
      figure <- as.character(localKeys[i,"FIGURE"])
      subarea <- as.character(localKeys[i,"SUBAREA"])
      area <- as.character(localKeys[i,"AREA"])
      # TO DO: remove subarea whole from probe features
      area_subarea <- paste(area,"_", ifelse (subarea=="","WHOLE",subarea) , sep="")
      source_pivot_filename <- pivot_file_name(marker, figure, "POSITION","")
      dest_pivot_filename <- pivot_file_name(marker, figure, area,subarea)

      # i <- 1
      if (!file.exists(dest_pivot_filename))
      {
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " File does not exists: ", dest_pivot_filename)
        {
          probe_features <- probe_features_get(area_subarea)
          probe_features$CHR <- as.factor(paste0("chr", probe_features$CHR))

          # annotate file
          pivot <- readr::read_delim(source_pivot_filename,
            col_types = readr::cols(
              .default = readr::col_double(),
              CHR = readr::col_character(),
              START = readr::col_integer(),
              END = readr::col_integer()
            ),
            show_col_types=FALSE, progress=FALSE)
          if(!is.null(pivot))
          {
            if(any(grepl("END", colnames(probe_features)))) #dplyr::inner_join
              pivot <- merge(probe_features,pivot, by = c("CHR", "START","END"), all.y = TRUE)
            else
              pivot <- merge(probe_features,pivot, by = c("CHR", "START"), all.y = TRUE)

            pivot <-subset(pivot, !is.na(eval(parse(text=area_subarea))))

            if(nrow(pivot)!=0)
            {
              pivot_colnames <- colnames(pivot)
              pivot_colnames[which(pivot_colnames==area_subarea)] <- "AREA"
              colnames(pivot) <- pivot_colnames

              pivot$CHR <- as.factor(pivot$CHR)

              probe_features <- probe_features[(probe_features$CHR %in% unique((pivot$CHR))), ]
              droplevels(probe_features$CHR)
              droplevels(pivot$CHR)

              colname_to_preserve <- !(colnames(pivot) %in%  c("PROBE","CHR","START","END","K27","K450","K850"))
              pivot <- pivot[, colname_to_preserve] %>% dplyr::distinct()
              # sum rows with the same AREA
              if(localKeys[i,"DISCRETE"])
                pivot <- pivot %>% dplyr::group_by(AREA) %>% dplyr::summarise_all(list(sum = sum))
              else
                pivot <- pivot %>% dplyr::group_by(AREA) %>% dplyr::summarise_all(list(mean = mean))
              readr::write_delim(pivot, dest_pivot_filename, col_names = T, progress = FALSE)
              rm(pivot)
              #
            }
          }
          if(ssEnv$showprogress)
            progress_bar(sprintf("Annotating position pivots."))
        }
      }
    }


  end_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotation genomic areas file finished in ", difftime(end_time,start_time,units = "mins")," minutes.")
  #
}

