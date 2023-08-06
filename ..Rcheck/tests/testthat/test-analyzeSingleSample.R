testthat::test_that("analyze_single_sample",{

  ssEnv <-  init_env(tempFolder, parallel_strategy = "sequential")
  ssEnv <- get_session_info()


  get_meth_tech(methylation_data)

  sp <- analyze_single_sample( values = methylation_data[,1],
                      sliding_window_size = sliding_window_size,
                      thresholds = thresholds,
                      figure = "HYPO",
                      sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")] ,
                      bonferroni_threshold = bonferroni_threshold,
                      probe_features = probe_features)


  outputFolder <- dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPO"))
  fileName <- file_path_build(outputFolder,c(Sample_ID,"MUTATIONS","HYPO"), "bed")
  testthat::expect_true(file.exists(fileName))

  unlink(tempFolder)
  close_env()

})
