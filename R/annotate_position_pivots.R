#' @title create an annotated file for each marker, figure, area and subarea, each file has all the sample_groups used to calculate epimutation
#'
#' @return nothing
#' @importFrom doRNG %dorng%
annotate_position_pivots <- function ()
{
  start_time <- Sys.time()
  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
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

  # doesn't work with parallel, tests throws error
  for(i in 1:nrow(localKeys))
    # foreach::foreach(i=1:nrow(localKeys), .export = variables_to_export) %dorng%
  {
    marker <- as.character(localKeys[i,"MARKER"])
    figure <- as.character(localKeys[i,"FIGURE"])
    subarea <- as.character(localKeys[i,"SUBAREA"])
    area <- as.character(localKeys[i,"AREA"])
    # TO DO: remove subarea whole from probe features
    area_subarea <- paste(area,"_", ifelse (subarea=="","WHOLE",subarea) , sep="")
    source_pivot_filename <- pivot_file_name_parquet(marker, figure, "POSITION","")
    dest_pivot_filename <- pivot_file_name_parquet(marker, figure, area,subarea)
    # i <- 1
    if (!file.exists(dest_pivot_filename))
    {
      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " File does not exists: ", dest_pivot_filename)
      {
        probe_features <- probe_features_get(area_subarea)
        probe_features$CHR <- as.factor(paste0("chr", probe_features$CHR))
        # probe_features <-subset(probe_features, !is.na(eval(parse(text=area_subarea))))

        # annotate file
        if(file.exists(source_pivot_filename))
        {
          log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotatinbg, reading pivot.")

          pivot <- polars::pl$read_parquet(source_pivot_filename)$to_data_frame()
          if(any(grepl("END", colnames(probe_features)))) #dplyr::inner_join
            pivot <- merge(probe_features,pivot, by = c("CHR", "START","END"), all.y = TRUE)
          else
            pivot <- merge(probe_features,pivot, by = c("CHR", "START"), all.y = TRUE)

          log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating, annotaion executed.")
          pivot <-subset(pivot, !is.na(eval(parse(text=area_subarea))))

          if(nrow(pivot)!=0)
          {
            pivot_colnames <- colnames(pivot)
            pivot_colnames[which(pivot_colnames==area_subarea)] <- "AREA"
            colnames(pivot) <- pivot_colnames

            # pivot$CHR <- as.factor(pivot$CHR)
            # probe_features <- probe_features[(probe_features$CHR %in% unique((pivot$CHR))), ]

            existing_cols <- names(pivot)
            cols_to_remove <- c("PROBE","CHR","START","END","K27","K450","K850")
            cols_to_remove <- cols_to_remove[cols_to_remove %in% existing_cols]

            pivot <- polars::as_polars_df(pivot)$drop(cols_to_remove)

            # droplevels(probe_features$CHR)
            # droplevels(pivot$CHR)

            # colname_to_preserve <- !(colnames(pivot) %in%  c("PROBE","CHR","START","END","K27","K450","K850"))
            # pivot <- pivot[, colname_to_preserve] %>% dplyr::distinct()

            # sum rows with the same AREA
            if (localKeys[i, "DISCRETE"]) {
              # pivot <- pivot %>%
              #   dplyr::group_by(AREA) %>%
              #   dplyr::summarise(across(everything(), sum, .names = "{.col}"))
              pivot <- pivot$group_by("AREA", maintain_order = TRUE)$sum()
            } else {
              # pivot <- pivot %>%
              #   dplyr::group_by(AREA) %>%
              #   dplyr::summarise(across(everything(), mean, .names = "{.col}"))
              pivot <- pivot$group_by("AREA", maintain_order = TRUE)$mean()
            }
            log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotatinbg, aggregation executed.")


            # pivot <- polars::as_polars_df(pivot)
            pivot$write_parquet(dest_pivot_filename)
            # arrow::write_parquet(pivot,dest_pivot_filename)
            log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotatinbg, writing executed marker %s figure %s, area %s subarea %s",marker, figure ,area, subarea)
            rm(pivot)
            #
          }
          else
            ssEnv$key_missed_areas_subareas <- unique(rbind(ssEnv$key_missed_areas_subareas, localKeys[i,c("AREA","SUBAREA")]))
        }
      }
    }

    if(ssEnv$showprogress)
      progress_bar(sprintf("Annotating position pivots."))
  }

  # remove missed keys
  selector <- !((ssEnv$keys_areas_subareas_markers_figures$AREA %in% ssEnv$key_missed_areas_subareas$AREA) & (ssEnv$keys_areas_subareas_markers_figures$SUBAREA %in% ssEnv$key_missed_areas_subareas$SUBAREA))
  ssEnv$keys_areas_subareas_markers_figures  <- ssEnv$keys_areas_subareas_markers_figures[selector,]

  update_session_info(ssEnv)
  end_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotation genomic areas file finished in ", difftime(end_time,start_time,units = "mins")," minutes.")
  gc()
  #
}

