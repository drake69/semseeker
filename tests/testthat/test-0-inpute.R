test_that("inpute", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  ex_ante <- sum(is.na(signal_data))
  nrow_ex_ante <- nrow(signal_data)

  ####################################################################################

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute="median")
  signal_data_median <- semseeker:::inpute_missing_values(signal_data)

  nrows_ex_post <- nrow(signal_data_median)
  ex_post <- sum(is.na(signal_data_median))

  testthat::expect_true(ex_ante!=0)
  testthat::expect_true(ex_post==0)
  testthat::expect_true(nrow_ex_ante - nrows_ex_post == nrow_missed)

  ####################################################################################

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute="mean")
  signal_data_mean <- semseeker:::inpute_missing_values(signal_data)

  nrows_ex_post <- nrow(signal_data_mean)
  ex_post <- sum(is.na(signal_data_mean))

  testthat::expect_true(ex_ante!=0)
  testthat::expect_true(ex_post==0)
  testthat::expect_true(nrow_ex_ante - nrows_ex_post == nrow_missed)

  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

