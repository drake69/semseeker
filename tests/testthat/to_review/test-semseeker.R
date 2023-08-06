test_that("semeeker", {


  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker( sample_sheet =  mySampleSheet_batch,methylation_data =  methylation_data_batch, result_folder = tempFolder,
             parallel_strategy = parallel_strategy, areas="GENE", figures="BOTH")

  tempresult_folder <- file.path(tempFolder,"Data","Control","MUTATIONS_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(localFileRes_both)>0)

  ####################################################################################

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAS_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(localFileRes_both)>0)

  ####################################################################################

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)

  ####################################################################################

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTARQ_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTARQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)

  ####################################################################################
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

