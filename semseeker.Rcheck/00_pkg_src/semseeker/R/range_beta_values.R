#--- range_beta_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of beta values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution
#' @importFrom doRNG %dorng%

range_beta_values <- function(populationMatrix, iqrTimes = 3) {

  # populationMatrixDim <- dim(populationMatrix)
  populationMatrix <- populationMatrix[, !(colnames(populationMatrix) %in% "PROBE")]
  beta_values <- populationMatrix
  row.names(beta_values) <- rownames(populationMatrix)
  # message("INFO: ", Sys.time(), " Starting beta thresholds calculation.")
  # progress_bar <- progressr::progressor(along = 1:nrow(beta_values))

  export = c("progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","beta_values","iqrTimes")
  r <- 1
  # nrow(beta_values)
  # for(r in 1:1000)
  # result <- foreach::foreach(r = 1:nrow(beta_values), .combine = "rbind", .export = export) %dorng%
  #   {
  #     beta_row <- beta_values[r,]
  #     b_values <- as.vector(t(beta_row))
  #     temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
  #     betaQ1Values <-  temp[1]
  #     betaQ3Values <- temp[2]
  #     beta_median_values <- temp[3]
  #     # b_values <- as.vector(t(beta_values[r,]))
  #     # betaQ1Values <-  stats::quantile(b_values, 0.25)
  #     # betaQ3Values <- stats::quantile(b_values, 0.75)
  #     # beta_median_values <- stats::quantile(b_values, 0.5)
  #     betaValuesIQR <- stats::IQR(b_values)
  #
  #     beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  #     beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))
  #
  #     temp_result <- data.frame("beta_inferior_thresholds"= beta_inferior_thresholds,
  #                               "beta_superior_thresholds"= beta_superior_thresholds,
  #                               "beta_median_values"= beta_median_values)
  #     row.names(temp_result) <- row.names(b_values)
  #     # colnames(temp_result) <- c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values")
  #     # if(attr(r, "nb_cycles") %% 10 == 0)
  #     #   progress_bar(sprintf("probe#: %s",names(b_values)))
  #     temp_result
  #   }

  message("INFO: ", Sys.time(), " Starting beta thresholds calculation.")
  # progress_bar_2 <- progressr::progressor(along = 1:nrow(beta_values)/10)
  # external_var <- 0

  # define a function to increase the external variable
  # increase_var <- function(x){
  #   external_var <<- external_var + x
  #   if(external_var %% 20 == 0)
  #     progress_bar_2()
  # }
  th <- future.apply::future_apply(beta_values,1,  get_th <- function(beta_row)
  {
    # message(names(beta_row))
    # message(colnames(beta_row))
    # message(rownames(beta_row))
    b_values <- as.vector(t(beta_row))
    temp <- stats::quantile(b_values, c(0.25,0.5,0.75))
    betaQ1Values <-  temp[1]
    betaQ3Values <- temp[2]
    beta_median_values <- temp[3]
    # betaQ1Values <-  stats::quantile(b_values, 0.25)
    # betaQ3Values <- stats::quantile(b_values, 0.75)
    # beta_median_values <- stats::quantile(b_values, 0.5)
    betaValuesIQR <- stats::IQR(b_values)

    beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
    # beta_inferior_thresholds[beta_inferior_thresholds<0] <- 0

    beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

    temp_result <- c("beta_inferior_thresholds"= beta_inferior_thresholds,
      "beta_superior_thresholds"= beta_superior_thresholds,
      "beta_median_values"= beta_median_values,
      "iqr" = betaValuesIQR,
      "q1" = betaQ1Values,
      "q3" = betaQ3Values)
    names(temp_result) <- names(b_values)
    # if(attr(beta_row, "nb_cycles") %% 10 == 0)
    #   progress_bar_2()
    # increase_var(1)
    # external_var <<- external_var + 1
    # if(external_var %% 10 == 0)
    #   progress_bar_2()
    temp_result
  })


  result <- as.data.frame(t(th))
  colnames(result) <- c("beta_inferior_thresholds","beta_superior_thresholds", "beta_median_values","iqr","q1","q3")
  # message("\n")
  # colnames(values) <- c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values")
  # row.names(values) <- row.names(beta_median_values)


  # result <- list(beta_inferior_thresholds = beta_inferior_thresholds,
  #                beta_superior_thresholds = beta_superior_thresholds,
  #                beta_median_values = beta_median_values)

  message("INFO: ", Sys.time(), " Thresholds defined for: ", nrow(result$beta_inferior_thresholds), " probe_features.")
  return(result)
}
