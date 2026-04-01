inference_file_name <- function(inference_detail, marker)
{
  ssEnv <- get_session_info()

  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
  family_test <- inference_detail$family_test
  transformation <- inference_detail$transformation
  independent_variable <- gsub(" ","", inference_detail$independent_variable)
  depth_analysis <- inference_detail$depth_analysis

  file_result_prefix <- paste(depth_analysis, as.character(independent_variable),sep = "_")
  if(is.null(covariates) || length(covariates)  ==  0)
    file_suffix <- ""
  else
    file_suffix <- paste(covariates, collapse = "_")

  fileNameResults <- file_path_build(ssEnv$result_folderInference,c( marker,file_result_prefix , as.character(transformation), as.character(family_test), file_suffix),"csv")

}
