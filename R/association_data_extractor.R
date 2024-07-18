association_data_extractor <- function(inference_details,destination_folder="", result_folder="", ...)
{

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  if(result_folder!="")
    ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- get_session_info()

  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  final_results <- data.frame()
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
        results_inference <- filter_sql(inference_detail$areas_sql_condition, results_inference)

        # remove any column with name containg SAMPLES_SQL_CONDITION
        results_inference <- results_inference[,!grepl("SAMPLES_SQL_CONDITION", colnames(results_inference))]
        utils::write.csv2(results_inference, fileNameResults, row.names = FALSE)

        if(inference_detail$samples_sql_condition!="")
          results_inference$SAMPLES_SQL_CONDITION <- inference_detail$samples_sql_condition

        final_results <- plyr::rbind.fill(final_results, results_inference)
        # log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," sql executed:", areas_sql_condition)
        filePathResults <- file.path(destination_folder, paste0(inference_detail$family_test, "_", markers[a], ".csv"))
        if(destination_folder != "")
          utils::write.csv2(results_inference, filePathResults, row.names = FALSE)
      }
    }
  }
  return(final_results)
}


summary_association_analysis <- function(inference_details,destination_folder="", result_folder="", ...)
{
  if(result_folder!="")
    ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- get_session_info()
  association_data <- association_data_extractor(inference_details, destination_folder, result_folder, ...)


  available_metrics <- toupper(semseeker::metrics_properties[,"Metric"])

  if(any(grepl("PVALUE",colnames(association_data))))
    available_metrics <- c(available_metrics, colnames(association_data)[grepl("PVALUE",colnames(association_data))])

  # remove not numeric columns metrics
  metrics_to_remove <- colnames(association_data[,!sapply(association_data, is.numeric)])
  available_metrics <- available_metrics[!(available_metrics %in% metrics_to_remove)]

  available_metrics <- available_metrics[available_metrics %in% colnames(association_data)]

  #  create a summary table for the association analysis grouping by AREA,SUBAREA,MARKER,FIGURE and SAMPLES_SQL_CONDITION if exists
  if(any("SAMPLES_SQL_CONDITION" %in% colnames(association_data))) {
    summary_table <- association_data %>%
      dplyr::group_by(AREA, SUBAREA, MARKER, FIGURE, SAMPLES_SQL_CONDITION) %>%
      dplyr::summarise(dplyr::across(available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  } else {
    summary_table <- association_data %>%
      dplyr::group_by(AREA, SUBAREA, MARKER, FIGURE) %>%
      dplyr::summarise(dplyr::across(available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  }



  summary_table <- as.data.frame(summary_table)
  return(summary_table)
}
