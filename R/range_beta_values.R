#--- range_beta_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of beta values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution
#' @importFrom doRNG %dorng%

range_beta_values <- function(populationMatrix, iqrTimes = 3) {

  ssEnv <- get_session_info()

  # populationMatrixDim <- dim(populationMatrix)
  populationMatrix <- populationMatrix[, !(colnames(populationMatrix) %in% "PROBE")]
  beta_values <- populationMatrix
  row.names(beta_values) <- rownames(populationMatrix)
  # message("INFO: ", Sys.time(), " Starting beta thresholds calculation.")
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(beta_values))

  message("INFO: ", Sys.time(), " Starting beta thresholds calculation.")
  export = c("progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","beta_values","iqrTimes","ssEnv")
  r <- 1

  if(ncol(beta_values) > 500)
  {
    # for(r in 1:1000)
    result <- foreach::foreach(r = 1:nrow(beta_values), .combine = "rbind", .export = export) %dorng%
      {
        beta_row <- beta_values[r,]
        b_values <- as.vector(t(beta_row))
        temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
        betaQ1Values <-  temp[1]
        betaQ3Values <- temp[2]
        beta_median_values <- temp[3]
        betaValuesIQR <- stats::IQR(b_values)

        beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
        beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

        temp_result <- data.frame("beta_inferior_thresholds"= beta_inferior_thresholds,
          "beta_superior_thresholds"= beta_superior_thresholds,
          "beta_median_values"= beta_median_values,
          "iqr" = betaValuesIQR,
          "q1"= betaQ1Values,
          "q3"= betaQ3Values)
        row.names(temp_result) <- row.names(b_values)
        # colnames(temp_result) <- c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values")
        # if(ssEnv$showprogress)
        #   progress_bar(sprintf("%s",names(b_values)))
        temp_result
      }

    gc()
  }
  else
  {
    th <- future.apply::future_apply(beta_values, 1 ,  get_th <- function(beta_row)
    {
      b_values <- as.vector(t(beta_row))
      temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
      betaQ1Values <-  temp[1]
      betaQ3Values <- temp[2]
      beta_median_values <- temp[3]
      betaValuesIQR <- stats::IQR(b_values)

      beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
      beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

      temp_result <- c("beta_inferior_thresholds"= beta_inferior_thresholds,
        "beta_superior_thresholds"= beta_superior_thresholds,
        "beta_median_values"= beta_median_values,
        "iqr" = betaValuesIQR,
        "q1" = betaQ1Values,
        "q3" = betaQ3Values)
      names(temp_result) <- names(b_values)
      temp_result
    })
    result <- as.data.frame(t(th))
  }

  colnames(result) <- c("beta_inferior_thresholds","beta_superior_thresholds", "beta_median_values","iqr","q1","q3")
  message("INFO: ", Sys.time(), " Thresholds defined for: ", nrow(result$beta_inferior_thresholds), " probe_features.")
  return(result)
}
