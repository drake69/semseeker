test_that("get-meth_tech", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")


  init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "BOTH", anomalies = "DELTAS", metaareas = "GENE")

  probes <- semseeker::PROBES

  ####################################################################################

  methylation_data_27 <- subset(methylation_data, rownames(methylation_data) %in% probes[probes$k27,"PROBE"])
  ssEnv <- get_meth_tech(methylation_data_27)
  testthat::expect_true(ssEnv$tech=="k27")

  ####################################################################################

  methylation_data_450 <- subset(methylation_data, rownames(methylation_data) %in% probes[probes$k450,"PROBE"])
  ssEnv <- get_meth_tech(methylation_data_450)
  testthat::expect_true(ssEnv$tech=="k450")

  ####################################################################################

  methylation_data_850 <- subset(methylation_data, rownames(methylation_data) %in% probes[probes$k850,"PROBE"])
  ssEnv <- get_meth_tech(methylation_data_850)
  testthat::expect_true(ssEnv$tech=="k850")

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  close_env()

})
