#' @export
covarage_analysis_report <- function (result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  signal_data_path <- pivot_file_name_parquet("SIGNAL","MEAN","PROBE","")
  if (!file.exists(signal_data_path))
  {
    log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " Signal data is missing")
    close_env()
    stop()
  }
  signal_data <- polars::pl$read_parquet(signal_data_path)$to_data_frame()
  ssEnv <- get_meth_tech(signal_data)
  if (ssEnv$tech=="WGBS")
  {
    log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " WGBS data is not supported for coverage analysis")
    close_env()
    stop()
  }
  coverage_analysis(signal_data$AREA)
  close_env()
}
