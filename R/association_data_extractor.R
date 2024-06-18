association_data_extractor <- function(inference_details, sql_conditions=c(),destination_folder, result_folder, ...)
{

  
  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  for(z in 1:nrow(inference_details))
  {
    start_time <- Sys.time()
    inference_detail <- inference_details[z,]
    for (a in 1:length(markers) )
    {
      if (exists("results_inference"))
        rm(results_inference)

      keys <- localKeys[localKeys$MARKER==markers[a],]
      keys <- unique(keys)
      fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference)
      if (file.exists(fileNameResults))
      {
        results_inference <- utils::read.csv2(fileNameResults, header  =  T)
        results_inference <- subset(results_inference, FIGURE %in% keys$FIGURE)
        results_inference <- filter_sql(sql_conditions, results_inference)

        log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," sql executed:", sql)
        filePathResults <- file.path(destination_folder, paste0(inference_detail$family_test, "_", markers[a], ".csv"))
        utils::write.csv2(results_inference, filePathResults, row.names = FALSE)
      }
    }
  }
}
