filter_metrics <- function(metrics, transformation){

  affected_by_tranformation <- toupper(c("MAE", "RMSE","MSLE", "SSE","R_SQUARED_ADJ","R_SQUARED","MSE",
    "MAE_TEST", "RMSE_TEST","MSLE_TEST", "SSE_TEST","R_SQUARED_ADJ_TEST","R_SQUARED_TEST","MSE_TEST",
    "STD.ERROR","EFFECT_SIZE_ESTIMATE","pinball_loss","POWER","STATISTIC_PARAMETER","kw_runk_sum"))

  #
  if (transformation!="scale")
  {
    to_remove_metrics <- metrics %in% affected_by_tranformation
    metrics <- metrics[!to_remove_metrics]
  }
  metrics <- sort(unique(c(metrics)))

  return(metrics)
}
