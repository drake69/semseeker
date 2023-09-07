testthat::test_that("delta_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  delta_single_sample(
    values = methylation_data[,1],
    high_thresholds = beta_superior_thresholds,
    low_thresholds = beta_inferior_thresholds,
    sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")],
    beta_medians = beta_medians,
    probe_features = probe_features
  )

  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
  outputFolder <- semseeker:::dir_check_and_create(result_folderData,c("Control","DELTAS_BOTH"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAS","BOTH"), "bedgraph")
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()

})
