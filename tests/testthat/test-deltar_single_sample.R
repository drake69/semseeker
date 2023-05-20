testthat::test_that("deltar_single_sample",{

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- semseeker:::init_env(tempFolder)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  semseeker:::deltar_single_sample(
    values = methylation_data[,1],
    high_thresholds = beta_superior_thresholds,
    low_thresholds = beta_inferior_thresholds,
    sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")],
    iqr = iqr,
    probe_features = probe_features
  )

  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
  outputFolder <- semseeker:::dir_check_and_create(result_folderData,c("Control","DELTAR_BOTH"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAR","BOTH"), "bedgraph")
  # deltar <- read.csv(fileName, sep="\t")
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()

})
