test_that("test_init_env", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <<- tempFolders[1]

  # test minimum number of parameters
  ssEnv <- semseeker:::init_env(tempFolder)
  testthat::expect_true(ssEnv$result_folderData == paste0(tempFolder,"/Data"))
  semseeker:::close_env()

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  testthat::expect_true(ssEnv$maxResources == "10")
  semseeker:::close_env()

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,areas = c("GENE","DMR"), showprogress=TRUE, start_fresh=TRUE)
  testthat::expect_true(length(unique(ssEnv$keys_areas_subareas[,"AREA"]))!=1)

  unlink(tempFolder,recursive = TRUE)
  assign("ssEnv", NULL, envir=.pkgglobalenv)

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "ERROR:")
  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, markers="HYPPO"), "ERROR:")
  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "ERROR:")
  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, subareas="HYPPO"), "ERROR:")
  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, test_it="HYPPO"), "ERROR:")

  ####################################################################################

}
)

test_that("inpute", {

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

  signal_data <<- signal_data_median
})

test_that("signal_range_values", {

  # test range calculation with missed values throws error
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, start_fresh=TRUE)
  testthat::expect_error( semseeker:::signal_range_values(signal_data), "^ERROR:")

  ####################################################################################
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute = "median", start_fresh=FALSE)
  signal_thresholds <<- semseeker:::signal_range_values(signal_data)

  ####################################################################################
  # check output has expected column names
  testthat::expect_true(sum(colnames(signal_thresholds) %in% c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values","iqr","q1","q3"))==6)
  # test no probe are lost
  testthat::expect_true(nrow(signal_thresholds)==nrow(signal_data))
  # test lower limits are different from upper limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_superior_thresholds[1])
  # test median is different from lower limits
  testthat::expect_true(signal_thresholds$signal_inferior_thresholds[1]!=signal_thresholds$signal_median_values[1])
  # test median is different from IQR
  testthat::expect_true(signal_thresholds$iqr[1]!=signal_thresholds$signal_median_values[1])
  # test median is different from Q1
  testthat::expect_true(signal_thresholds$q1[1]!=signal_thresholds$signal_median_values[1])
  # test median is different from Q3
  testthat::expect_true(signal_thresholds$q3[1]!=signal_thresholds$signal_median_values[1])
  # test median is different from upper limits
  testthat::expect_true(sum(signal_thresholds$signal_inferior_thresholds<signal_thresholds$signal_superior_thresholds)/nrow(signal_thresholds)>0.9)
  # test thresholds row are the same of probes of input matrix
  testthat::expect_true( nrow(signal_thresholds)==nprobes)
  # test median is lower than Q3
  testthat::expect_true(sum(signal_thresholds$signal_median_values<signal_thresholds$q3)/nrow(signal_thresholds)>0.9)
  # test median is higher than Q1
  testthat::expect_true(sum(signal_thresholds$signal_median_values>signal_thresholds$q1)/nrow(signal_thresholds)>0.9)

})



