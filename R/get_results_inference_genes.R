get_results_inference_genes <- function (inference_details,marker, pvalue = 0.05, adjust_per_area = F, adjust_globally = F,pvalue_column="PVALUEADJ_ALL_BH",adjustment_method = "BH")
{

  ssEnv <- get_session_info()
  resultFolder <- ssEnv$result_folderInference

  inferenceFile <- inference_file_name(inference_detail, marker)
  if(adjust_per_area && adjust_globally)
  {
    message("Can adjust per area or globbaly not both!")
    stop()
  }

  results_inference <- read.csv(paste(resultFolder,inferenceFile, sep=""), sep=";", dec=",")
  results_inference <- results_inference[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","BETA",pvalue_column)]
  results_inference <- subset(results_inference,AREA=="GENE" & AREA_OF_TEST!="TOTAL")
  results_inference_hh <- results_inference
  results_inference_hh$FIGURE <- "HYPER_HYPO"
  results_inference <- rbind(results_inference, results_inference_hh)
  if(adjust_globally)
    results_inference[,pvalue_column] <- p.adjust(results_inference$PVALUE, method = adjustment_method)

  seq <- 0

  # check if convariates are used
  if(length(unique(results_inference$COVARIATES))>1)
    applied_model <- paste(c(unique(results_inference$INDIPENDENT.VARIABLE), unlist(strsplit(unique(results_inference$COVARIATES)," "))), collapse = "_")
  else
    applied_model <- paste(c(unique(results_inference$INDIPENDENT.VARIABLE)), collapse = "_")

  transformation <- unique(results_inference$transformation)
  FAMILY <- unique(results_inference$FAMILY)
  adjustment_text = "no_adjustment"
  if(adjust_per_area)
    adjustment_text= "adjusted_per_area"
  if(adjust_globally)
    adjustment_text= "adjusted_globally"
  # markers <- unique(keys$MARKER)

  if(adjust_per_area)
    for (a in unique(results_inference$AREA))
    {
      results_inference_area <- results_inference[results_inference$AREA==a,]
      results_inference_area[,pvalue_column] <- p.adjust(results_inference_area$PVALUE, method = adjustment_method)
      results_inference[results_inference$AREA==a,] <- results_inference_area
    }

  results_inference <- subset(results_inference, results_inference[,pvalue_column] < pvalue)
  results_inference <- na.omit(results_inference)

  return(results_inference)
}
