metrics_ranking <- function (metric,data_frame, figure, column_to_rank ="REBASED"){

  #
  the_lower_the_better_markers <- toupper(c("COUNT_MISSED","MAE", "RMSE","MSLE", "SSE","MAPE","MPE","MSE",
    "MAE_TEST", "RMSE_TEST","MSLE_TEST", "SSE_TEST","MAPE_TEST","MPE_TEST","MSE_TEST","STD.ERROR","BELOW_QUANTILE",
    "PINBALL_LOSS","QQ_DISTANCE"))

  # R_SQUARED	R_SQUARED_ADJ
  the_higher_the_better_markers <- toupper(c("R_SQUARED_ADJ","R_SQUARED","COUNT_SIGN","Edata_frameECT_SIZE_ESTIMATE","Edata_frameECT_SIZE_MAGNITUDE",
    "RBC", "POWER", "STATISTIC_PARAMETER","JSD","QQ_CORRELATION",
    "RANK_BESERIAL_CORRELATION","Wilcox_Value","kw_runk_sum","R_SQUARED_ADJ_TEST","R_SQUARED_TEST"))


  if(grepl("PVALUE",metric))
    the_lower_the_better_markers <- c(the_lower_the_better_markers, metric)

  data_frame[,column_to_rank] <- abs(data_frame[,column_to_rank])
  if(metric %in% the_lower_the_better_markers){
    data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = T),]
  }else if(metric %in% the_higher_the_better_markers){
    data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = F),]
  }else{
    log_event("ERROR: Metric not found")
    stop()
  }
  # create a new column with the rank
  data_frame$RANK <- 1:nrow(data_frame)

  # get only duplicates of column_to_rank
  duplicates <- data_frame[duplicated(data_frame[,column_to_rank]) | duplicated(data_frame[,column_to_rank], fromLast = T),column_to_rank]
  duplicates <- unique(duplicates)

  for (r in 1:length(duplicates)){
    duplicate_value <- duplicates[r]
    # assign the same rank where REBASED is the same
    data_frame[data_frame[,column_to_rank] == duplicate_value,"RANK"] <- data_frame[r,"RANK"]
  }

  # # add the marker and figure to the data frame
  # scores <- rbind(scores, )
  scores <- data.frame("MARKER"=data_frame$MARKER,"FIGURE"=figure,"METRIC"=metric,"SCORE"=data_frame$RANK)

  return(scores)
}
