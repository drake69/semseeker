  test_that("analyze_population", {

  init_env(tempFolder, parallel_strategy = "sequential")

  sp <- analyze_population(methylation_data=methylation_data,
                    sliding_window_size = 11,
                    sliding_window_size = sliding_window_size,
                    beta_thresholds = beta_thresholds,
                    sample_sheet = mySampleSheet,
                    bonferroni_threshold = bonferroni_threshold,
                    probe_features = probe_features,
                    bonferroni_threshold = 0.01,
                    )

  sp$Sample_Group <- mySampleSheet$Sample_Group
  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))

  unlink(tempFolder,recursive = TRUE)
  close_env()

})

