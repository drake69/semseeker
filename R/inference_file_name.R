inference_file_name <- function(inference_detail, marker, folder,file_extension="csv", prefix = "", suffix="")
{
  ssEnv <- get_session_info()

  #
  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
  covariates_dummy <- inference_detail$covariates_dummy
  covariates_dummy <- if(length(covariates_dummy) !=  0 && !is.null(covariates_dummy)) unlist(t(strsplit( gsub(" ","",covariates_dummy),split  =  "+", fixed  =  T)))
  covariates_pca <- ifelse(is.null(inference_detail$covariates_pca),FALSE,inference_detail$covariates_pca)

  family_test <- inference_detail$family_test
  transformation <- inference_detail$transformation
  independent_variable <- gsub(" ","", inference_detail$independent_variable)
  depth_analysis <- inference_detail$depth_analysis

  file_suffix <- ""
  file_result_prefix <- paste(as.character(independent_variable),as.character(transformation),as.character(family_test),sep="_")
  if(!is.null(covariates) & (length(covariates) !=0))
  {
    long_covariates <- length(covariates) > 2
    # split each covariates by _
    if (long_covariates)
    {
      covariates <- unlist(t(strsplit( gsub(" ","",covariates),split  =  "_", fixed  =  T)))
      covariates <- unique(covariates)
    }
    file_suffix <- paste(covariates, collapse = "_")
  }

  if (!is.null(covariates_dummy) && (length(covariates_dummy)  !=  0))
  {
    long_covariates <- length(covariates_dummy) > 2
    # split each covariates_dummy by _
    if (long_covariates)
    {
      covariates_dummy <- unlist(t(strsplit( gsub(" ","",covariates_dummy),split  =  "_", fixed  =  T)))
      covariates_dummy <- unique(covariates_dummy)
    }
    file_suffix <- c(file_suffix, paste(covariates_dummy, collapse = "_"))
  }

  if(length(covariates_dummy))
    file_suffix <- c(file_suffix, paste("dummy", sep = "_"))

  if(covariates_pca)
    file_suffix <- c(file_suffix, paste("pca", sep = "_"))

  file_result_prefix <- paste("DEPTH",depth_analysis, file_result_prefix,sep = "_")
  file_result_prefix = paste(prefix, file_result_prefix, sep = "_")
  file_suffix = paste(file_suffix, suffix, sep = "_")

  if(length(inference_detail$samples_sql_condition)>0)
    folder <- dir_check_and_create(folder,name_cleaning(inference_detail$samples_sql_condition))


  fileNameResults <- file_path_build(folder,c(as.character(marker),file_result_prefix ,file_suffix) , file_extension)

  return(fileNameResults)
}


