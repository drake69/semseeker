inference_file_name <- function(inference_detail, marker, folder,file_extension="csv", prefix = "", suffix="")
{
  ssEnv <- semseeker:::get_session_info()

  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
  family_test <- inference_detail$family_test
  transformation <- inference_detail$transformation
  independent_variable <- gsub(" ","", inference_detail$independent_variable)
  depth_analysis <- inference_detail$depth_analysis

  file_result_prefix <- paste(as.character(independent_variable),as.character(transformation),as.character(family_test),sep="_")
  if(is.null(covariates) || length(covariates)  ==  0)
    file_suffix <- ""
  else
    file_suffix <- paste(covariates, collapse = "_")


  file_result_prefix <- paste("DEPTH",depth_analysis, file_result_prefix,sep = "_")
  file_result_prefix = paste(prefix, file_result_prefix, sep = "")
  file_suffix = paste(file_suffix, suffix, sep = "")

  fileNameResults <- semseeker:::file_path_build(folder,c(as.character(marker),file_result_prefix ,file_suffix) , file_extension)

  return(fileNameResults)
}


