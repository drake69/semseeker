test_that("range_beta_values", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################
  beta_thresholds_local <- semseeker:::range_beta_values(methylation_data, iqrTimes = iqrTimes)
  testthat::expect_true(sum(colnames(beta_thresholds_local)==c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values","iqr","q1","q3"))==6)

  ####################################################################################
  testthat::expect_true(nrow(beta_thresholds_local)==nrow(methylation_data))
  ####################################################################################
  testthat::expect_true(beta_thresholds_local$beta_inferior_thresholds[1]!=beta_thresholds_local$beta_superior_thresholds[1])
  ####################################################################################
  testthat::expect_true(beta_thresholds_local$beta_inferior_thresholds[1]!=beta_thresholds_local$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(beta_thresholds_local$iqr[1]!=beta_thresholds_local$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(beta_thresholds_local$q1[1]!=beta_thresholds_local$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(beta_thresholds_local$q3[1]!=beta_thresholds_local$beta_median_values[1])
  ####################################################################################
  testthat::expect_true(all(beta_thresholds_local$beta_inferior_thresholds<beta_thresholds_local$beta_superior_thresholds))

  ####################################################################################
  testthat::expect_true( nrow(beta_thresholds_local)==nprobes)


  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

