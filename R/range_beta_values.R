#--- range_beta_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of beta values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation
#'
#' @param iqrTimes inter quartile ratio used to normalize
#'
#' @return methylation matrix as normalized distribution
#' @importFrom doRNG %dorng%

range_beta_values <- function(populationMatrix, iqrTimes = 3) {

  populationMatrixDim <- dim(populationMatrix)
  beta_values <- as.data.frame(populationMatrix[, 2:populationMatrixDim[2]])
  row.names(beta_values) <- populationMatrix[, 1]

  r <- 1
  # nrow(beta_values)
  # for(r in 1:1000)
  result <- foreach::foreach(r = 1:nrow(beta_values), .combine = "rbind", .export = c("beta_values","iqrTimes")) %dorng%
    {
      b_values <- as.vector(t(beta_values[r,]))
      betaQ1Values <-  stats::quantile(b_values, 0.25)
      betaQ3Values <- stats::quantile(b_values, 0.75)
      beta_median_values <- stats::quantile(b_values, 0.5)
      betaValuesIQR <- stats::IQR(b_values)

      beta_inferior_thresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
      beta_superior_thresholds <- (betaQ3Values + (iqrTimes * betaValuesIQR))

      temp_result <- data.frame("beta_inferior_thresholds"= beta_inferior_thresholds,
                                "beta_superior_thresholds"= beta_superior_thresholds,
                                "beta_median_values"= beta_median_values)
      row.names(temp_result) <- row.names(b_values)
      # colnames(temp_result) <- c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values")
      temp_result
    }

  # colnames(values) <- c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values")
  # row.names(values) <- row.names(beta_median_values)


  # result <- list(beta_inferior_thresholds = beta_inferior_thresholds,
  #                beta_superior_thresholds = beta_superior_thresholds,
  #                beta_median_values = beta_median_values)

  message("Thresholds defined for: ", nrow(result$beta_inferior_thresholds), " probes.")
  return(result)
}
