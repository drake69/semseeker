metrics_ranking <- function (metric,data_frame, figure, column_to_rank ="REBASED"){

  #
  the_lower_the_better_markers <- toupper(semseeker::metrics_properties[semseeker::metrics_properties$Higher_the_Better==FALSE,"Metric"])

  # R_SQUARED	R_SQUARED_ADJ
  the_higher_the_better_markers <-  toupper(semseeker::metrics_properties[semseeker::metrics_properties$Higher_the_Better==TRUE,"Metric"])


  if(grepl("PVALUE",metric))
    the_lower_the_better_markers <- c(the_lower_the_better_markers, metric)

  metric <- toupper(metric)
  # data_frame[,column_to_rank] <- abs(data_frame[,column_to_rank])
  if(metric %in% the_lower_the_better_markers){
    data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = T),]
  }else if(metric %in% the_higher_the_better_markers){
    data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = F),]
  }else{
    log_event("ERROR: Metric not found")
    stop()
  }

  # browser()
  values_to_rank <- unique(data_frame[,column_to_rank])
  values_to_rank <- data.frame("VALUE"=values_to_rank)
  # create a new column with the rank
  # values_to_rank$RANK <- 1:nrow(values_to_rank)

  values_to_rank$RANK <- cut(values_to_rank$VALUE, breaks = 100, right = FALSE, labels = FALSE)
  # merge the rank with the data frame
  data_frame <- merge(data_frame, values_to_rank, by.x=column_to_rank, by.y="VALUE", all.x=TRUE)

  # # add the marker and figure to the data frame
  # scores <- rbind(scores, )
  scores <- data.frame("MARKER"=data_frame$MARKER,"FIGURE"=figure,"METRIC"=metric,"SCORE"=data_frame$RANK)

  return(scores)
}
