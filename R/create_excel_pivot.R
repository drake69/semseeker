#' @importFrom doRNG %dorng%
create_excel_pivot <-  function() {

  start_time <- Sys.time()
  ssEnv <- get_session_info()
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  # # count the number of files created previously
  # existing_files <- foreach::foreach(k = 1:nrow(localKeys),.combine = sum) %dorng%
  #   {
  #     area <-  as.character(localKeys[k,"AREA"])
  #     subarea <-  as.character(localKeys[k,"SUBAREA"])
  #     marker <- as.character(localKeys[k,"MARKER"])
  #     figure <-  as.character(localKeys[k,"FIGURE"])
  #     pivot_file_name <-  as.character(localKeys[k,"COMBINED"])
  #     pivot_subfolder <- dir_check_and_create(reportFolder, marker)
  #     pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_file_name,extension =  ".csv" ,add_gz=TRUE)
  #     if (file.exists(pivot_file_name)){
  #       1
  #     } else {
  #       0
  #     }
  #   }
  #
  #
  #
  # # exists only the signal pivot saved at the beginning of the elaboration
  # if(existing_files == 1)
  #   ssEnv$keys_areas_subareas_markers_figures_missed <- NULL
  #
  # if (!is.null(ssEnv$keys_areas_subareas_markers_figures_missed))
  #   # remove the missed keys from the localKeys
  #   localKeys <- localKeys[!(localKeys$COMBINED %in% ssEnv$keys_areas_subareas_markers_figures_missed$COMBINED),]

  # sample_sheet <- utils::read.csv2(file.path(ssEnv$result_folderData,"sample_sheet_result.csv"))
  # create_multiple_bed(sample_sheet)
  # annotate_bed()

  if(ssEnv$showprogress)
    progress_bar_ann <- progressr::progressor(along=1:nrow(localKeys))
  else
    progress_bar_ann <- ""

  variables_to_export <- c("ssEnv", "reportFolder", "dir_check_and_create","localKeys",
    "progress_bar_ann","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","file_path_build", "read_annotated_bed")

  log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," Creating Excel ", nrow(localKeys) ,  " Pivots")
  for (k in 1:nrow(localKeys))
  # created <- foreach::foreach(k = 1:nrow(localKeys), .export = variables_to_export, .combine = max) %dorng%
    {
      # k <- 400
      area <-  as.character(localKeys[k,"AREA"])
      subarea <-  as.character(localKeys[k,"SUBAREA"])
      marker <- as.character(localKeys[k,"MARKER"])
      figure <-  as.character(localKeys[k,"FIGURE"])
      pivot_file_name <-  as.character(localKeys[k,"COMBINED"])

      pivot_subfolder <- dir_check_and_create(reportFolder, marker)
      pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_file_name,extension =  ".csv" ,add_gz=TRUE)
      if (!file.exists(pivot_file_name))
      {
        annotatedData <-  read_annotated_bed(figure,marker,area,subarea)
        annotatedData <- subset(annotatedData, annotatedData$VALUE != 0 )

        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y") , " Creating for marker:", marker ," and figure:", figure)

        if(!plyr::empty(annotatedData))
        {
          if(marker=="SIGNAL")
          {
            annotatedData <- reshape2::dcast(data = annotatedData, formula = SAMPLEID + SAMPLE_GROUP ~ AREA, value.var = "VALUE",
              fun.aggregate = mean, drop = TRUE)
          }
          else
            annotatedData <- reshape2::dcast(data = annotatedData, formula =  SAMPLEID + SAMPLE_GROUP ~ AREA, value.var = "VALUE",
              fun.aggregate = sum, drop = TRUE)

          utils::write.table(t(annotatedData), gzfile(pivot_file_name), row.names = T, col.names = F, sep=";")
          end_time <- Sys.time()
          log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y") , " Finished to create for marker:",
            marker ," figure:", figure, " area: ", area, " marker ", marker,
            " - Time taken: ", difftime(end_time,start_time, units = "mins") , " minutes.")
        }
      }
      if(ssEnv$showprogress)
        progress_bar_ann(sprintf("Creating pivot file."))
      # k
    }
    created <-nrow(localKeys)

  end_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y") , " Finished to create ", created,  " excel pivot - Time taken: ", difftime(end_time,start_time, units = "mins") , " minutes.")
}
