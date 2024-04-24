#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sample_sheet dataframe with at least a column Sample_ID to identify samples
#' @param signal_data matrix of methylation data
#' @param result_folder where the result will be saved
#' @param ... other options to execute elaborations
#'
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#' @importFrom doRNG %dorng%
#'
semseeker <- function(sample_sheet,
                      signal_data,
                      result_folder,
                      ... ) {

  init_env( result_folder= result_folder, ...)
  ssEnv <- get_session_info()
  write.csv2(sample_sheet, file = file.path(ssEnv$result_folderData, "sample_sheet_original.csv"), row.names = FALSE)

  # set digits to 22
  withr::local_options(list(digits = 22))
  sliding_window_size <- 11

  if(is.data.frame(sample_sheet) & is.data.frame(signal_data))
  {
    sample_sheet <-list(sample_sheet)
    signal_data <- list(signal_data)
  } else
  {
    if(is.data.frame(sample_sheet) | is.data.frame(signal_data))
      stop("both sample_sheet and signal_data should be data frame!")
    if(length(sample_sheet)!=length(signal_data))
      stop("both sample_sheet and signal_data should have been list with the same length!")
  }

  if(length(signal_data)>1)
  {
    d <- 1
    for(d in 1:length(signal_data))
    {
      if (ssEnv$signal_intrasample)
        probes_to_preserve <- row.names((signal_data[[d]]))
      else
        probes_to_preserve <- row.names(stats::na.omit(signal_data[[d]]))
    }
  }
  else
    probes_to_preserve <- row.names(signal_data[[1]])

  batch_id <- 1
  for(batch_id in 1:length(sample_sheet))
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Working on batch:",batch_id)
    sample_sheet_local <- sample_sheet[[batch_id]]
    signal_intrasample <- TRUE
    if (ssEnv$signal_intrasample)
      signal_data_local <- signal_data[[batch_id]]
    else
      signal_data_local <- stats::na.omit(signal_data[[batch_id]])
    signal_data_local <- signal_data_local[rownames(signal_data_local) %in% probes_to_preserve,]
    sample_sheet_local <- analyze_batch(signal_data_local, sample_sheet_local, batch_id)
    if(exists("sample_sheet_result"))
      sample_sheet_result <- plyr::rbind.fill(sample_sheet_result, sample_sheet_local)
    else
      sample_sheet_result <- sample_sheet_local
    utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  }

  sample_sheet <- sample_sheet_result
  utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Saving Sample Sheet with Results! ")

  if(length(sample_sheet$Sample_Group=="Reference")>0)
    sample_groups <- c("Reference","Control","Case")
  else
    sample_groups <- c("Control","Case")

  # seems to change the ssEnv reducing the items TODO
  ssEnv <- get_session_info()
  ssEnv$keys_sample_groups <- data.frame("SAMPLE_GROUP"=sample_groups)
  update_session_info(ssEnv)

  annotate_bed()
  create_excel_pivot()

 log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  "Starting inference Analysis.")
  # inferenceAnalysis(ssEnv$result_folderData = ssEnv$result_folderData, ssEnv$session_folder= ssEnv$session_folder, inferenceDetails)
  # future::autoStopCluster(computationCluster)
  # doFuture::stopImplicitCluster()

  # geneontology_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData, fileName = fileName)
  # euristic_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData)

  # if(length(signal_data)>1)
  #   batch_correlation_check(ssEnv)

  close_env()
}
