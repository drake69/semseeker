#' @importFrom doRNG %dorng%
create_excel_pivot <-  function() {

  start_time <- Sys.time()
  ssEnv <- get_session_info()
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  localKeys <- ssEnv$keys_areas_subareas_markers_figures
  # areas <- ssEnv$keys_areas

  sample_sheet <- utils::read.csv2(file.path(ssEnv$result_folderData,"sample_sheet_result.csv"))
  # sample_sheet <- deltaq_get(sample_sheet)
  # sample_sheet <- deltarq_get(sample_sheet)
  # write.csv2(sample_sheet, file.path(ssEnv$result_folderData,"sample_sheet_result.csv"), row.names = FALSE)

  create_multiple_bed(sample_sheet)
  annotate_bed()

  if(ssEnv$showprogress)
    progress_bar_ann <- progressr::progressor(along=1:nrow(localKeys))
  else
    progress_bar_ann <- ""

  variables_to_export <- c("ssEnv", "annotatedData", "subGroupLabel", "SAMPLE_GROUP", "reportFolder", "area",
    "sheetList","dir_check_and_create","localKeys",
    "progress_bar_ann","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","file_path_build", "read_annotated_bed")

  for (k in 1:length(localKeys))
  # foreach::foreach(k = 1:length(localKeys), .export = variables_to_export) %dorng%
  {
    area <-  as.character(localKeys[k,"AREA"])
    subarea <-  as.character(localKeys[k,"SUBAREA"])
    marker <- as.character(localKeys[k,"MARKER"])
    figure <-  as.character(localKeys[k,"FIGURE"])
    pivot_file_name <- localKeys$COMBINED[k]

    pivot_subfolder <- dir_check_and_create(reportFolder, marker)
    fileName <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_file_name,extension =  ".csv" ,add_gz=TRUE)
    if (!file.exists(fileName))
    {
      annotatedData <-  read_annotated_bed(figure,marker,area,subarea)
      annotatedData <- subset(annotatedData, annotatedData$VALUE != 0 )

      # browser()
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

          utils::write.table(t(annotatedData), gzfile(fileName), row.names = T, col.names = F, sep=";")
      }
    }
    if(ssEnv$showprogress)
      progress_bar_ann(sprintf("Creating pivot file."))
  }
  gc()
  end_time <- Sys.time()
  log_event(paste0("INFO: ", end_time, " Finished to create excel pivot - Time taken: ", end_time - start_time))
}
