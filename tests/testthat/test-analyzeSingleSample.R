testthat::test_that("analyze_single_sample",{

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <-  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)
  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  sp <- semseeker:::analyze_single_sample( values = methylation_data[,1],
                      sliding_window_size = sliding_window_size,
                      thresholds = rep(1,nrow(methylation_data)),
                      figure = "HYPO",
                      sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")] ,
                      bonferroni_threshold = 0.05,
                      probe_features = probe_features)


  outputFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPO"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"MUTATIONS","HYPO"), "bed")
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
