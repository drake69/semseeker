test_that("get-meth_tech", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "BOTH", markers = "DELTAS", areas = "GENE")

  probe_features <- semseeker::pp_tot

  ####################################################################################

  methylation_data_27 <- subset(methylation_data, rownames(methylation_data) %in% probe_features[probe_features$K27,"PROBE"])
  ssEnv <- semseeker:::get_meth_tech(methylation_data_27)
  testthat::expect_true(ssEnv$tech=="K27")

  ####################################################################################

  methylation_data_450 <- subset(methylation_data, rownames(methylation_data) %in% probe_features[probe_features$K450,"PROBE"])
  ssEnv <- semseeker:::get_meth_tech(methylation_data_450)
  testthat::expect_true(ssEnv$tech=="K450")

  ####################################################################################

  methylation_data_850 <- subset(methylation_data, rownames(methylation_data) %in% probe_features[probe_features$K850,"PROBE"])
  ssEnv <- semseeker:::get_meth_tech(methylation_data_850)
  testthat::expect_true(ssEnv$tech=="K850")

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()

})
