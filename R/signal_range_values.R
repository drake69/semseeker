#--- signal_range_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of signal values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @return methylation matrix as normalized distribution
#' @importFrom doRNG %dorng%

signal_range_values <- function(populationMatrix) {

  ssEnv <- get_session_info()

  # populationMatrixDim <- dim(populationMatrix)
  populationMatrix <- populationMatrix[, !(colnames(populationMatrix) %in% "PROBE")]
  signal_values <- populationMatrix
  min_values <- apply(signal_values, 1, min, na.rm=TRUE)
  max_values <- apply(signal_values, 1, max, na.rm=TRUE)
  row.names(signal_values) <- rownames(populationMatrix)
  log_event("DEBUG: ", Sys.time(), " Starting signal thresholds calculation.")
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(signal_values))

  log_event("INFO: ", Sys.time(), " Starting signal thresholds calculation.")
  export = c("progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_values","ssEnv")
  r <- 1

  if(ncol(signal_values) > 500)
  {
    # for(r in 1:1000)
    result <- foreach::foreach(r = 1:nrow(signal_values), .combine = "rbind", .export = export) %dorng%
      {
        signal_row <- signal_values[r,]
        signal_values <- as.vector(t(signal_row))
        temp <- stats::quantile(signal_values, c(0.25,0.5,0.75), na.rm=TRUE)
        signalQ1Values <-  temp[1]
        signalQ3Values <- temp[2]
        signal_median_values <- temp[3]
        signalValuesIQR <- stats::IQR(signal_values)

        signal_inferior_thresholds <- max((signalQ1Values - (ssEnv$iqrTimes * signalValuesIQR)), min(signal_values))
        signal_superior_thresholds <- min((signalQ3Values + (ssEnv$iqrTimes * signalValuesIQR)), max(signal_values))

        temp_result <- data.frame(
          "signal_inferior_thresholds"= signal_inferior_thresholds,
          "signal_superior_thresholds"= signal_superior_thresholds,
          "signal_median_values"= signal_median_values,
          "iqr" = signalValuesIQR,
          "q1"= signalQ1Values,
          "q3"= signalQ3Values)
        row.names(temp_result) <- row.names(signal_values)
        # temp_result$PROBE <- row.names(signal_values)
        # colnames(temp_result) <- c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values")
        if(ssEnv$showprogress)
          progress_bar(sprintf("%s",names(signal_values)))
        temp_result
      }

    gc()
  }
  else
  {
    th <- future.apply::future_apply(signal_values, 1 ,  get_th <- function(signal_row)
    {
      signal_values <- as.vector(t(signal_row))
      temp <- stats::quantile(signal_values, c(0.25,0.5,0.75), na.rm = TRUE)
      signalQ1Values <-  temp[1]
      signalQ3Values <- temp[2]
      signal_median_values <- temp[3]
      signalValuesIQR <- stats::IQR(signal_values)

      signal_inferior_thresholds <- max((signalQ1Values - (ssEnv$iqrTimes * signalValuesIQR)), min(signal_values))
      signal_superior_thresholds <- min((signalQ3Values + (ssEnv$iqrTimes * signalValuesIQR)), max(signal_values))

      temp_result <- c(
        "signal_inferior_thresholds"= signal_inferior_thresholds,
        "signal_superior_thresholds"= signal_superior_thresholds,
        "signal_median_values"= signal_median_values,
        "iqr" = signalValuesIQR,
        "q1" = signalQ1Values,
        "q3" = signalQ3Values)
      names(temp_result) <- names(signal_values)
      # temp_result$PROBE <- names(signal_values)
      temp_result
    })
    result <- as.data.frame(t(th))
  }

  colnames(result) <- c("signal_inferior_thresholds","signal_superior_thresholds", "signal_median_values","iqr","q1","q3")
  result$PROBE <- row.names(signal_values)
  # if (min_values==result$signal_inferior_thresholds)
  # {
  #   browser()
  # }
  # if (max_values==result$signal_superior_thresholds)
  # {
  #   browser()
  # }
  log_event("INFO: ", Sys.time(), " Thresholds defined for: ", nrow(result), " probe_features.")
  return(result)
}
