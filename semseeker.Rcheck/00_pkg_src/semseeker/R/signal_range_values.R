#--- signal_range_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of signal values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution
#' @importFrom doRNG %dorng%

signal_range_values <- function(populationMatrix, iqrTimes = 3) {

  # populationMatrixDim <- dim(populationMatrix)
  populationMatrix <- populationMatrix[, !(colnames(populationMatrix) %in% "PROBE")]
  signal_values <- populationMatrix
  row.names(signal_values) <- rownames(populationMatrix)
  # message("INFO: ", Sys.time(), " Starting signal thresholds calculation.")
  # progress_bar <- progressr::progressor(along = 1:nrow(signal_values))

  export = c("progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_values","iqrTimes")
  r <- 1
  # nrow(signal_values)
  # for(r in 1:1000)
  # result <- foreach::foreach(r = 1:nrow(signal_values), .combine = "rbind", .export = export) %dorng%
  #   {
  #     signal_row <- signal_values[r,]
  #     b_values <- as.vector(t(signal_row))
  #     temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
  #     signalQ1Values <-  temp[1]
  #     signalQ3Values <- temp[2]
  #     signal_median_values <- temp[3]
  #     # b_values <- as.vector(t(signal_values[r,]))
  #     # signalQ1Values <-  stats::quantile(b_values, 0.25)
  #     # signalQ3Values <- stats::quantile(b_values, 0.75)
  #     # signal_median_values <- stats::quantile(b_values, 0.5)
  #     signalValuesIQR <- stats::IQR(b_values)
  #
  #     signal_inferior_thresholds <- (signalQ1Values - (iqrTimes * signalValuesIQR))
  #     signal_superior_thresholds <- (signalQ3Values + (iqrTimes * signalValuesIQR))
  #
  #     temp_result <- data.frame("signal_inferior_thresholds"= signal_inferior_thresholds,
  #                               "signal_superior_thresholds"= signal_superior_thresholds,
  #                               "signal_median_values"= signal_median_values)
  #     row.names(temp_result) <- row.names(b_values)
  #     # colnames(temp_result) <- c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values")
  #     # if(attr(r, "nb_cycles") %% 10 == 0)
  #     #   progress_bar(sprintf("probe#: %s",names(b_values)))
  #     temp_result
  #   }

  message("INFO: ", Sys.time(), " Starting signal thresholds calculation.")
  # progress_bar_2 <- progressr::progressor(along = 1:nrow(signal_values)/10)
  # external_var <- 0

  # define a function to increase the external variable
  # increase_var <- function(x){
  #   external_var <<- external_var + x
  #   if(external_var %% 20 == 0)
  #     progress_bar_2()
  # }
  th <- future.apply::future_apply(signal_values,1,  get_th <- function(signal_row)
  {
    # message(names(signal_row))
    # message(colnames(signal_row))
    # message(rownames(signal_row))
    b_values <- as.vector(t(signal_row))
    temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
    signalQ1Values <-  temp[1]
    signalQ3Values <- temp[2]
    signal_median_values <- temp[3]
    # signalQ1Values <-  stats::quantile(b_values, 0.25)
    # signalQ3Values <- stats::quantile(b_values, 0.75)
    # signal_median_values <- stats::quantile(b_values, 0.5)
    signalValuesIQR <- stats::IQR(b_values)

    signal_inferior_thresholds <- (signalQ1Values - (iqrTimes * signalValuesIQR))
    # signal_inferior_thresholds[signal_inferior_thresholds<0] <- 0

    signal_superior_thresholds <- (signalQ3Values + (iqrTimes * signalValuesIQR))

    temp_result <- c("signal_inferior_thresholds"= signal_inferior_thresholds,
      "signal_superior_thresholds"= signal_superior_thresholds,
      "signal_median_values"= signal_median_values,
      "iqr" = signalValuesIQR,
      "q1" = signalQ1Values,
      "q3" = signalQ3Values)
    names(temp_result) <- names(b_values)
    # if(attr(signal_row, "nb_cycles") %% 10 == 0)
    #   progress_bar_2()
    # increase_var(1)
    # external_var <<- external_var + 1
    # if(external_var %% 10 == 0)
    #   progress_bar_2()
    temp_result
  })


  result <- as.data.frame(t(th))
  colnames(result) <- c("signal_inferior_thresholds","signal_superior_thresholds", "signal_median_values","iqr","q1","q3")
  # message("\n")
  # colnames(values) <- c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values")
  # row.names(values) <- row.names(signal_median_values)


  # result <- list(signal_inferior_thresholds = signal_inferior_thresholds,
  #                signal_superior_thresholds = signal_superior_thresholds,
  #                signal_median_values = signal_median_values)

  message("INFO: ", Sys.time(), " Thresholds defined for: ", nrow(result$signal_inferior_thresholds), " probe_features.")
  return(result)
}
