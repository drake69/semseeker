test_that("range_beta_values", {


  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################
  rr <- range_beta_values(methylation_data, iqrTimes = iqrTimes)
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

  unlink(tempFolder,recursive = TRUE)
  close_env()
})

