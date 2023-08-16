#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sample_sheet dataframe with at least a column Sample_ID to identify samples
#' @param methylation_data matrix of methylation data
#' @param bonferroni_threshold = 0.05 #threshold to define which pValue adjusted to define an epilesion
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param iqrTimes how many times below the first quartile and over the third quartile the interqauartile is "added" to define the outlier
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param result_folder where the result will be saved
#' @param ... other options to filter elaborations
#'
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#' @importFrom doRNG %dorng%
#'
semseeker <- function(sample_sheet,
                      methylation_data,
                      result_folder,
                      bonferroni_threshold = 0.05,
                      maxResources = 90,
                      iqrTimes = 3 ,
                      parallel_strategy ="multisession",
                      ... ) {


  unlink(result_folder, recursive = TRUE)
  init_env( result_folder= result_folder, maxResources= maxResources, parallel_strategy = parallel_strategy, ...)
  ssEnv <- get_session_info()

  # set digits to 22
  withr::local_options(list(digits = 22))
  sliding_window_size <- 11

  if(is.data.frame(sample_sheet) & is.data.frame(methylation_data))
  {
    sample_sheet <-list(sample_sheet)
    methylation_data <- list(methylation_data)
  } else
  {
    if(is.data.frame(sample_sheet) | is.data.frame(methylation_data))
      stop("both sample_sheet and methylation_data should be data frame!")
    if(length(sample_sheet)!=length(methylation_data))
      stop("both sample_sheet and methylation_data should have been list with the same length!")
  }

  if(length(methylation_data)>1)
  {
    d <- 1
    for(d in 1:length(methylation_data))
    {
      if(exists("probes_to_preserve"))
        probes_to_preserve <- probes_to_preserve[ probes_to_preserve %in% row.names(stats::na.omit(methylation_data[[d]]))]
      else
        probes_to_preserve <- row.names(stats::na.omit(methylation_data[[d]]))
    }
  }
  else
    probes_to_preserve <- row.names(methylation_data[[1]])

  batch_id <- 1
  for(batch_id in 1:length(sample_sheet))
  {
    message("INFO: ", Sys.time(), " Working on batch:",batch_id)
    sample_sheet_local <- sample_sheet[[batch_id]]
    methylation_data_local <- stats::na.omit(methylation_data[[batch_id]])
    methylation_data_local <- methylation_data_local[rownames(methylation_data_local)%in%probes_to_preserve,]
    sample_sheet_local <- analyze_batch(methylation_data_local, sample_sheet_local, sliding_window_size, bonferroni_threshold,iqrTimes, batch_id)
    if(exists("sample_sheet_result"))
      sample_sheet_result <- plyr::rbind.fill(sample_sheet_result, sample_sheet_local)
    else
      sample_sheet_result <- sample_sheet_local
    utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  }

  sample_sheet <- sample_sheet_result
  utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  message("INFO: ", Sys.time(), " Saving Sample Sheet with Results! ", Sys.time())

  if(length(sample_sheet$Sample_Group=="Reference")>0)
    sample_groups <- c("Reference","Control","Case")
  else
    sample_groups <- c("Control","Case")

  ssEnv$keys_sample_groups <- sample_groups
  update_session_info(ssEnv)

  annotate_bed()
  create_excel_pivot()

  # message("Starting inference Analysis.")
  # inferenceAnalysis(ssEnv$result_folderData = ssEnv$result_folderData, ssEnv$session_folder= ssEnv$session_folder, inferenceDetails)
  # future::autoStopCluster(computationCluster)
  # doFuture::stopImplicitCluster()

  # geneontology_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData, fileName = fileName)
  # euristic_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData)

  # if(length(methylation_data)>1)
  #   batch_correlation_check(ssEnv)

  close_env()
}
