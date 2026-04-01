# get result data of an association analyisis and repeat over significant areas
#' @export
association_to_association <- function(inference_details_origin, inference_details,result_folder,
  maxResources = 90, parallel_strategy  = "multicore",start_fresh = FALSE, ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  for (a in seq_along(markers))
  {

    marker <- markers[a]

    inference_filename <- inference_file_name(inference_detail = inference_details, folder = ssEnv$result_folderInference, marker = marker, file_extension = "csv")

    if(file.exists(inference_filename))
      next

    inference_source <- association_results_get(inference_details_origin, marker, adjust_per_area = FALSE, adjust_globally = FALSE,
      pvalue_column="PVALUE_ADJ_ALL_FDR",adjustment_method = "BH", area ="GENE",
      omit_na = TRUE, significance = TRUE)

    areas <- as.vector(unique(inference_source$AREA_OF_TEST))

    association_analysis( inference_details = inference_details, result_folder = result_folder, areas_selection=areas,
      maxResources = maxResources,  parallel_strategy  = parallel_strategy,start_fresh = start_fresh,
      areas = unique(ssEnv$keys_areas_subareas_markers_figures$AREA), subareas = unique(ssEnv$keys_areas_subareas_markers_figures$SUBAREA),
      markers = marker, verbosity= ssEnv$verbosity, figures=unique(ssEnv$keys_areas_subareas_markers_figures$FIGURE),
      showprogress = ssEnv$showprogress)



    inference_filename_origin <- inference_file_name(inference_detail = inference_details_origin, folder = ssEnv$result_folderInference, marker = marker, file_extension = "csv")
    data_origin <- utils::read.csv2(inference_filename_origin)
    data <- utils::read.csv2(inference_filename)

    # remove from data origin where figure, area, subare and are_of_test are not in data
    data_origin$KEY <- paste(data_origin$FIGURE, data_origin$AREA, data_origin$SUBAREA, data_origin$AREA_OF_TEST, sep = "_")
    data$KEY <- paste(data$FIGURE, data$AREA, data$SUBAREA, data$AREA_OF_TEST, sep = "_")
    data_origin <- data_origin[!(data_origin$KEY %in% data$KEY),]
    data_origin$FAMILY_TEST <- inference_details$family_test

    data <- rbind(data, data_origin)
    data$KEY <- NA
    association_analysis_save_results(data, inference_filename, inference_details$family_test, FALSE)

  }


}
