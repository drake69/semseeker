#' @export
recover <- function(result_folder, maxResources = 90,
  parallel_strategy  = "multicore",start_fresh = FALSE, ...)
{

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)

  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will trye to recover the analysys for project \n in ", ssEnv$result_folderData)
  summary_file <- file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  if(!file.exists(summary_file))
    summary_file <- file_path_build( ssEnv$result_folderData, "1_SAMPLE_SHEET_ORIGINALS","csv")
  if(!file.exists(summary_file))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No sample sheet found in ", ssEnv$result_folderData)
    return()
  }
  study_summary <-   utils::read.csv2(summary_file)
  study_summary <- create_position_pivots(study_summary,ssEnv$keys_markers_figures)
  study_summary <- deltaX_get(study_summary)
  summary_file <- file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  utils::write.csv2(study_summary,summary_file, row.names = FALSE)
  annotate_position_pivots()

}
