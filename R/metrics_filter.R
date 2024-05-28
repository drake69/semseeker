metrics_filter <- function(metrics, transformation){

  affected_by_tranformation <- toupper(metrics_properties[metrics_properties$Affected_by_Scaling==TRUE,"Metric"])

  if (transformation!="scale")
  {
    to_remove_metrics <- metrics %in% affected_by_tranformation
    metrics <- metrics[!to_remove_metrics]
  }
  metrics <- sort(unique(c(metrics)))

  return(metrics)

}
