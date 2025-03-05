#' Calculate Stochastic Epigenetic Mutations from a Methylation Dataset
#'
#' The `semseeker` function processes a methylation dataset to identify stochastic epigenetic mutations and generates output reports. This involves working with sample sheets and signal data to produce pivot tables and bedgraph files.
#'
#' @param sample_sheet A data frame containing at least a column named `Sample_ID` to identify samples. This can be a single data frame or a list of data frames.
#' @param signal_data A matrix of methylation data. This can be a single data frame or a list of data frames.
#' @param result_folder A string specifying the directory where the results will be saved.
#' @param ... Additional arguments for further processing options, including:
#'   - `parallel_strategy`: Strategy for parallel execution. Possible values are `none`, `multisession`, `sequential`, `multicore`, and `cluster`.
#'   - `maxResources`: Percentage of available cores to be used, default is 90 percent, rounded to the lowest integer.
#'   - `signal_intrasample`: A logical value indicating whether the signal data is intrasample. Default is `FALSE`.
#'   - `sliding_window_size`: An integer specifying the size of the sliding window. Default is 11.
#'   - `alpha`: A numeric value specifying the alpha threshold for the analysis. Default is 0.05.
#'   - `showprogress`: A logical value indicating whether to show a progress bar. Default is `TRUE`.
#'   - `iqrTimes`: A numeric value specifying the interquartile range multiplier to identify aberration in the data. Default is 3.
#'   - `sex_chromosome_remove`: A logical value indicating whether to remove the sex chromosomes. Default is `TRUE`.
#'   - `plot_format`: A string specifying the plot format. Default is "png".
#'   - `verbosity`: A numeric value specifying the verbosity level. Default is 0.
#'   - `marker`: A vector of string specifying the marker to be used for the analysis. Default are "MUTATIONS", "LESIONS","DELTARQ","DELTAQ","DELTAS","DELTAR","MEAN","DELTARP","DELTAP"
#'   - `areas`: A vector of string specifying the areas to be used for the analysis. Default are "GENE","CHR","ISLAND","PROBE"
#'
#' @return The function writes multiple files to the specified `result_folder`, including the processed sample sheet and result files in CSV format, along with pivot tables and bedgraph files.
#' @export
#' @importFrom doRNG %dorng%
#'
semseeker <- function(sample_sheet,
  signal_data,
  result_folder,
  ... ) {

  init_env( result_folder= result_folder, ...)

  ssEnv <- get_session_info()
  log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will search MARKERS for project \n in ", ssEnv$result_folderData)

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
    # browser()
    sample_sheet_local <- sample_sheet[[batch_id]]
    sample_sheet_local$Sample_ID <- name_cleaning(sample_sheet_local$Sample_ID)
    utils::write.csv2(sample_sheet_local, file = file_path_build(ssEnv$result_folderData, paste0(batch_id,"_sample_sheet_original"),"csv",FALSE))
    signal_intrasample <- TRUE
    signal_data_local <- signal_data[[batch_id]]
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " working on batch:", batch_id, " of ", nrow(signal_data_local), " rows and ", ncol(signal_data_local), " samples.")
    colnames(signal_data_local) <- name_cleaning(colnames(signal_data_local))
    signal_data_local <- signal_data_local[rownames(signal_data_local) %in% probes_to_preserve,]
    signal_data_local <- substitute_infinite(signal_data_local)
    gc()
    signal_data_local <- inpute_missing_values(signal_data_local)
    gc()
    signal_data_local <- stats::na.omit(signal_data_local)
    analyze_batch(signal_data_local, sample_sheet_local, batch_id)
    create_position_pivots(sample_sheet_local,ssEnv$keys_markers_figures)
    study_summary_total()
  }

  annotate_position_pivots()
  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " Saving Sample Sheet with Results! ")

  close_env()
}
