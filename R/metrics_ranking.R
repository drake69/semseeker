metrics_ranking <- function (metric,data_frame, scores, figure){

  # browser()
  the_lower_the_better_markers <- c("COUNT_MISSED","MAE", "RMSE","MSLE", "SSE","R_SQUARED_ADJ","R_SQUARED","MAPE","MPE","MSE",
    "MAE_TEST", "RMSE_TEST","MSLE_TEST", "SSE_TEST","R_SQUARED_ADJ_TEST","R_SQUARED_TEST","MAPE_TEST","MPE_TEST",
    "R_SQUARED_COMPL_TEST","R_SQUARED_ADJ_COMPL_TEST","MSE_TEST",
    "STD.ERROR")

  # R_SQUARED	R_SQUARED_COMPL	R_SQUARED_ADJ	R_SQUARED_ADJ_COMPL
  the_higher_the_better_markers <- c("COUNT_SIGN","EFFECT_SIZE_ESTIMATE","EFFECT_SIZE_MAGNITUDE", "RBC", "POWER", "STATISTIC_PARAMETER","JSD",
    "RANK_BESERIAL_CORRELATION","R_SQUARED_COMPL","R_SQUARED_ADJ_COMPL")


  if(grepl("PVALUE",metric))
    the_lower_the_better_markers <- c(the_lower_the_better_markers, metric)

  # browser()
  ff <- data_frame
  ff$REBASED <- abs(ff$REBASED)
  if(metric %in% the_lower_the_better_markers){
    ff <- ff[order(ff$REBASED, decreasing = T),]
  }else if(metric %in% the_higher_the_better_markers){
    ff <- ff[order(ff$REBASED, decreasing = F),]
  }else{
    stop("Metric not found")
  }
  # create a new column with the rank
  ff$RANK <- 1:nrow(ff)

  for (r in 1:nrow(ff)){
    # assign the same rank where REBASED is the same
    ff[ff$REBASED == ff[r,"REBASED"],"RANK"] <- ff[r,"RANK"]
  }

  # add the marker and figure to the data frame
  scores <- rbind(scores, data.frame("MARKER"=ff$MARKER,"FIGURE"=figure,"METRIC"=metric,"SCORE"=ff$RANK))

  return(scores)
}
