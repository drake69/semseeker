metrics_ranking <- function (metric,data_frame, figure, column_to_rank ="REBASED"){

  #
  the_lower_the_better_markers <- toupper(c("COUNT_MISSED","MAE", "RMSE","MSLE", "SSE","MAPE","MPE","MSE",
    "MAE_TEST", "RMSE_TEST","MSLE_TEST", "SSE_TEST","MAPE_TEST","MPE_TEST","MSE_TEST","STD.ERROR","BELOW_QUANTILE",
    "PINBALL_LOSS","QQ_DISTANCE"))

  # R_SQUARED	R_SQUARED_ADJ
  the_higher_the_better_markers <- toupper(c("R_SQUARED_ADJ","R_SQUARED","COUNT_SIGN","EFFECT_SIZE_ESTIMATE","EFFECT_SIZE_MAGNITUDE",
    "RBC", "POWER", "STATISTIC_PARAMETER","JSD","QQ_CORRELATION",
    "RANK_BESERIAL_CORRELATION","Wilcox_Value","kw_runk_sum","R_SQUARED_ADJ_TEST","R_SQUARED_TEST"))


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
  values_to_rank$RANK <- 1:nrow(values_to_rank)
  # merge the rank with the data frame
  data_frame <- merge(data_frame, values_to_rank, by.x=column_to_rank, by.y="VALUE", all.x=TRUE)

  # # add the marker and figure to the data frame
  # scores <- rbind(scores, )
  scores <- data.frame("MARKER"=data_frame$MARKER,"FIGURE"=figure,"METRIC"=metric,"SCORE"=data_frame$RANK)

  return(scores)
}
