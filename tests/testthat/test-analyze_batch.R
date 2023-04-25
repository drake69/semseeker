test_that("analize_batch", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  get_meth_tech(methylation_data)

  ####################################################################################

  batch_id <- 2
  sp <- analyze_batch(methylation_data =  methylation_data,
    sample_sheet =  mySampleSheet,
    sliding_window_size = sliding_window_size,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    batch_id = batch_id)


  sp$Sample_Group <- mySampleSheet$Sample_Group
  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))
  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_BOTH"])>0)>0)

  ####################################################################################
  unlink(tempFolder,recursive = TRUE)
  close_env()
})

