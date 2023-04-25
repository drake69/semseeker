  test_that("analize_population", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder, parallel_strategy = "sequential")

  sp <- analize_population(methylation_data=methylation_data,
                    sliding_window_size = sliding_window_size,
                    beta_superior_thresholds = beta_superior_thresholds,
                    beta_inferior_thresholds = beta_inferior_thresholds,
                    sample_sheet = mySampleSheet,
                    beta_medians = beta_medians,
                    bonferroni_threshold = bonferroni_threshold,
                    probe_features = probe_features
                    )

  sp$Sample_Group <- mySampleSheet$Sample_Group
  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))

  unlink(tempFolder,recursive = TRUE)
  close_env()

})

