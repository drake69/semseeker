test_that("analize_batch", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  batch_id <- 2
  sp <- semseeker:::analyze_batch(methylation_data =  methylation_data,
    sample_sheet =  mySampleSheet,
    sliding_window_size = sliding_window_size,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    batch_id = batch_id
  )


  sp$Sample_Group <- mySampleSheet$Sample_Group
  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))
  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_BOTH"])>0)>0)

  ####################################################################################
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

