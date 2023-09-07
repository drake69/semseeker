test_that("create_multiple_bed", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

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

  semseeker:::create_multiple_bed(mySampleSheet)
  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
  tempresult_folder <- semseeker:::dir_check_and_create(result_folderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(localFileRes_both)>0)

  ####################################################################################

  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()
})
