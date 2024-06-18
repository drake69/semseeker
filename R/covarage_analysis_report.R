#' @export
covarage_analysis_report <- function (signal_data, result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  signal_data <- as.data.frame(signal_data)
  get_meth_tech(signal_data)
  coverage_analysis(signal_data = signal_data)
  close_env()
}
