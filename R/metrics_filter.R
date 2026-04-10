metrics_filter <- function(metrics, transformation_y){

  if(transformation_y=="none")
    return(metrics)

  affected_by_tranformation <- toupper(SEMseeker::metrics_properties[SEMseeker::metrics_properties$Affected_by_Scaling==TRUE,"Metric"])

  if (transformation_y!="scale")
  {
    to_remove_metrics <- metrics %in% affected_by_tranformation
    metrics <- metrics[!to_remove_metrics]
  }
  metrics <- sort(unique(c(metrics)))

  return(metrics)

}
