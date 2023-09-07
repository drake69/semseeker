test_that("create_excel_pivot", {


  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  sp <- semseeker:::analyze_population(methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features

  )

  semseeker:::create_multiple_bed(sample_sheet = mySampleSheet)
  semseeker:::annotate_bed()

  ####################################################################################

  # create and read
  semseeker:::create_excel_pivot ()
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/GENE.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/PROBE.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/CHR.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/ISLAND.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/DMR.xlsx")))


  # unlink(tempFolder, recursive=TRUE)
  semseeker:::close_env()
})
