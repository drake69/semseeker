test_that("analize_batch", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder,
    parallel_strategy = parallel_strategy,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    inpute="median"
  )

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  batch_id <- 2
  sp <- semseeker:::analyze_batch(signal_data =  signal_data,
    sample_sheet =  mySampleSheet,
    batch_id = batch_id
  )


  sp$Sample_Group <- mySampleSheet$Sample_Group
  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))
  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_BOTH"])>0)>0)

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

