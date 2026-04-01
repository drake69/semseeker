test_that("analyze_population", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################


  
  sp <- semseeker:::analyze_population(signal_data=signal_data,
    sliding_window_size = sliding_window_size,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features
  )

  sp$Sample_Group <- mySampleSheet$Sample_Group

  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))


  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()

})

