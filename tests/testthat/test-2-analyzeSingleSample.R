test_that("analyze_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <-  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)
  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  sp <- semseeker:::analyze_single_sample(
                      values = signal_data[,1],
                      thresholds = signal_inferior_thresholds,
                      figure = "HYPO",
                      sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")] ,
                      probe_features = probe_features)


  outputFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPO"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"MUTATIONS","HYPO"), "bed", add_gz = TRUE)
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})
