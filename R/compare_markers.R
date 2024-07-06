compare_markers <- function(inference_detail){

  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  color_palette <- ssEnv$color_palette
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  aggregated_markers_results <- data.frame()
  for (a in 1:nrow(localKeys))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each study in studies
    for (s in 1:nrow(studies))
    {
      # get the inference details for the study
      temp_res <- get_results_areas_inference(inference_details = inference_detail, marker = MARKER, area= AREA,
         adjust_per_area = adjust_per_area, adjust_globally = adjust_globally, pvalue_column= pvalue_column,adjustment_method = adjustment_method, areas_sql_condition = "")
      if(nrow(temp_res) != 0)
        temp_res$STUDY <- studies[s,"STUDY"]
      aggregated_markers_results <- plyr::rbind.fill(aggregated_markers_results, temp_res)
    }
  }
}
