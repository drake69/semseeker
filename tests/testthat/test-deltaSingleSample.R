testthat::test_that("delta_single_sample",{

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- init_env(tempFolder)

  ####################################################################################

  get_meth_tech(methylation_data)

  ####################################################################################

  delta_single_sample(
    values = methylation_data[,1],
    high_thresholds = beta_superior_thresholds,
    low_thresholds = beta_inferior_thresholds,
    sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")],
    beta_medians = beta_medians,
    probe_features = probe_features
  )

  result_folderData  <-  dir_check_and_create(tempFolder, "Data")
  outputFolder <- dir_check_and_create(result_folderData,c("Control","DELTAS_BOTH"))
  fileName <- file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAS","BOTH"), "bedgraph")
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  close_env()

})
