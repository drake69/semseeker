test_that("signal_range_values", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################
  rr <- semseeker:::signal_range_values(signal_data, iqrTimes = iqrTimes)
  testthat::expect_true(sum(colnames(rr)==c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values"))==3)

  ####################################################################################
  testthat::expect_true(nrow(rr)==nrow(signal_data))
  ####################################################################################
  testthat::expect_true(rr$signal_inferior_thresholds[1]!=rr$signal_superior_thresholds[1])
  ####################################################################################
  testthat::expect_true(rr$signal_inferior_thresholds[1]!=rr$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$iqr[1]!=rr$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$q1[1]!=rr$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(rr$q3[1]!=rr$signal_median_values[1])
  ####################################################################################

  testthat::expect_true(all(rr$signal_inferior_thresholds<rr$signal_superior_thresholds))

  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

