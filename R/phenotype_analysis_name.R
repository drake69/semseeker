phenotype_analysis_name <- function(inference_detail, key,prefix = "", suffix="", pvalue_column, alpha, significance)
{
  ssEnv <- get_session_info()

  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
  family_test <- inference_detail$family_test
  transformation <- inference_detail$transformation
  independent_variable <- gsub(" ","", inference_detail$independent_variable)
  significance_label <- ifelse(significance, "significant", "non_significant")

  if(is.null(covariates) || length(covariates)  ==  0)
    file_suffix <- ""
  else
    file_suffix <- paste(covariates, collapse = "_")

  file_suffix = paste(file_suffix, suffix, sep = "")

  analysis_name <- paste(as.character(key$MARKER),as.character(key$FIGURE),as.character(key$AREA),as.character(key$SUBAREA), pvalue_column, significance_label, alpha, prefix , as.character(transformation), as.character(family_test), file_suffix, sep="_")

  return(analysis_name)
}
