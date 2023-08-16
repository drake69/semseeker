test_that("analyze_population", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################


  # browser()
  sp <- semseeker:::analyze_population(methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
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

