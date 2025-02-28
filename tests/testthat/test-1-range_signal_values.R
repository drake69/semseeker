test_that("signal_range_values", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]


  # test range calculation with missed values
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes)
  testthat::expect_error( semseeker:::signal_range_values(signal_data, batch_id), "^ERROR:")

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "median")
  signal_data <- semseeker:::inpute_missing_values(signal_data)

  ####################################################################################
  signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  testthat::expect_true(sum(colnames(signal_thresholds) %in% c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values","iqr","q1","q3"))==6)

  ####################################################################################
  # check thresholds file exists
  testthat::expect_true(file.exists(file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"parquet")))

  ####################################################################################
  # test no probe are lost
  testthat::expect_true(nrow(signal_thresholds)==nrow(signal_data))
  ####################################################################################
  # test lower limits are different from upper limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_superior_thresholds[1])
  ####################################################################################
  # test median is different from lower limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from IQR
  testthat::expect_true(signal_thresholds$iqr[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from Q1
  testthat::expect_true(signal_thresholds$q1[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from Q3
  testthat::expect_true(signal_thresholds$q3[1]!=signal_thresholds$signal_median_values[1])
  ####################################################################################
  # test median is different from upper limits
  testthat::expect_true(sum(signal_thresholds$signal_inferior_thresholds<signal_thresholds$signal_superior_thresholds)/nrow(signal_thresholds)>0.9)

  ####################################################################################
  # test thresholds row are the same of probes of input matrix
  testthat::expect_true( nrow(signal_thresholds)==nprobes)

  ####################################################################################
  # test median is lower than Q3
  testthat::expect_true(sum(signal_thresholds$signal_median_values<signal_thresholds$q3)/nrow(signal_thresholds)>0.9)

  ####################################################################################
  # test median is higher than Q1
  testthat::expect_true(sum(signal_thresholds$signal_median_values>signal_thresholds$q1)/nrow(signal_thresholds)>0.9)

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

