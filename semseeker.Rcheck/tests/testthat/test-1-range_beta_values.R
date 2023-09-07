test_that("range_beta_values", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################
  rr <- semseeker:::range_beta_values(methylation_data, iqrTimes = iqrTimes)
  testthat::expect_true(sum(colnames(rr)==c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values"))==3)

  ####################################################################################
  testthat::expect_true(nrow(rr)==nrow(methylation_data))
  ####################################################################################
  testthat::expect_true(rr$beta_inferior_thresholds[1]!=rr$beta_superior_thresholds[1])
  ####################################################################################
  testthat::expect_true(rr$beta_inferior_thresholds[1]!=rr$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$iqr[1]!=rr$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$q1[1]!=rr$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$q3[1]!=rr$beta_median_values[1])
  ####################################################################################

  testthat::expect_true(all(rr$beta_inferior_thresholds<rr$beta_superior_thresholds))

  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

