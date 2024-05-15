get_results_areas_inference <- function (inference_details, marker, adjust_per_area = F, adjust_globally = F,
  pvalue_column="PVALUE_ADJ_ALL_BH",adjustment_method = "BH", area ="GENE",
  omit_na = TRUE, sql_condition = "", significance = TRUE)
{
  ssEnv <- get_session_info()
  resultFolder <- ssEnv$result_folderInference

  inferenceFile <- inference_file_name(inference_details, marker, ssEnv$result_folderInference)
  if(adjust_per_area && adjust_globally)
  {
    log_event("ERROR: Can adjust per area or globbaly not both!", format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  if(!file.exists(inferenceFile))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),  " Inference file does not exist: ", inferenceFile)
    return(data.frame())
  }


  results_inference <- read.csv2(inferenceFile, sep=";", dec=",", row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  colnames(results_inference) <- gsub("[.]", "_", colnames(results_inference))

  # check columns exist
  if((!pvalue_column %in% colnames(results_inference)))
  {
    log_event("DEBUG:", pvalue_column, " column does not exist in inference file: ", inferenceFile)
    return(data.frame())
  }

  metrics_name_collect(results_inference)

  results_inference <- subset(results_inference,DEPTH==3)
  # results_inference <- results_inference[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE","DEPTH")]

  results_inference[results_inference$AREA==area,"AREA_OF_TEST"] <- gsub(results_inference[results_inference$AREA==area,"AREA_OF_TEST"] , pattern="_", replacement="-")

  if(adjust_globally)
    results_inference[,pvalue_column] <- p.adjust(results_inference[, pvalue_column], method = adjustment_method)

  adjustment_text = "no_adjustment"
  if(adjust_per_area)
    adjustment_text= "adjusted_per_area"
  if(adjust_globally)
    adjustment_text= "adjusted_globally"
  # markers <- unique(keys$MARKER)

  areas <- unique(results_inference$AREA)
  if(adjust_per_area)
    for (a in 1:length(areas))
    {
      area <- areas[a]
      results_inference_area <- results_inference[results_inference$AREA==area,]
      results_inference_area[,pvalue_column] <- p.adjust(results_inference_area[, pvalue_column], method = adjustment_method)
      results_inference[results_inference$AREA==area,] <- results_inference_area
    }

  #
  # preserve only subareas selected
  results_inference <- results_inference[results_inference$SUBAREA %in% unique(ssEnv$keys_areas_subareas$SUBAREA),]
  if (significance)
    results_inference <- subset(results_inference, results_inference[,pvalue_column] < as.numeric(ssEnv$alpha))
  else
    results_inference <- subset(results_inference, results_inference[,pvalue_column] >= as.numeric(ssEnv$alpha))
  # if(omit_na)
  #   results_inference <- na.omit(results_inference)

  #
  results_inference <- filter_sql(sql_condition, results_inference)
  log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," inference file loaded:", inferenceFile)
  return(results_inference)
}
