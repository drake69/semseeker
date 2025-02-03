#' @export
marker_distribution_info <- function(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = FALSE, ...)

  annotate_bed()
  create_excel_pivot()

  marker_quantization_metric()
  marker_fit_to_theoretical_distribution()

  close_env()

}
