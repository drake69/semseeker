test_that("create_excel_pivot", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, showprogress = TRUE)

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  multiple <- semseeker:::create_multiple_bed(mySampleSheet)
  dq <- deltaq_get(mySampleSheet)
  drq <- deltarq_get(mySampleSheet)
  semseeker:::annotate_bed()

  ####################################################################################

  # create and read
  semseeker:::create_excel_pivot()
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/GENE.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/PROBE.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/CHR.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/ISLAND.xlsx")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/DMR.xlsx")))


  # unlink(tempFolder, recursive=TRUE)
  semseeker:::close_env()
})
