test_that("signal_range_values", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes)

  ####################################################################################
  signal_thresholds_local <- semseeker:::signal_range_values(signal_data)
  testthat::expect_true(sum(colnames(signal_thresholds_local) %in% c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values","iqr","q1","q3"))==6)

  ####################################################################################
  testthat::expect_true(nrow(signal_thresholds_local)==nrow(signal_data))
  ####################################################################################
  testthat::expect_true(signal_thresholds_local$signal_inferior_thresholds[1]!=signal_thresholds_local$signal_superior_thresholds[1])
  ####################################################################################
  testthat::expect_true(signal_thresholds_local$signal_inferior_thresholds[1]!=signal_thresholds_local$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(signal_thresholds_local$iqr[1]!=signal_thresholds_local$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(signal_thresholds_local$q1[1]!=signal_thresholds_local$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(signal_thresholds_local$q3[1]!=signal_thresholds_local$signal_median_values[1])
  ####################################################################################
  testthat::expect_true(all(signal_thresholds_local$signal_inferior_thresholds<signal_thresholds_local$signal_superior_thresholds))

  ####################################################################################
  testthat::expect_true( nrow(signal_thresholds_local)==nprobes)

  ####################################################################################
  testthat::expect_true(all(signal_thresholds_local$signal_median_values<signal_thresholds_local$q3))

  ####################################################################################
  testthat::expect_true(all(signal_thresholds_local$signal_median_values>signal_thresholds_local$q1))

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

