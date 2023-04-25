test_that("semeeker", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- init_env(tempFolder, parallel_strategy = "sequential")

  ###############################################################################################
  ###############################################################################################

  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,
             parallel_strategy = "sequential", anomalies="DELTAQ", metaareas="GENE", figures="BOTH")

  # batch_correlation_check(ssEnv)
  tempresult_folder <- file.path(tempFolder,"Data","Control","MUTATIONS_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  # localFileRes <- read.table(fileToRead, sep="\t")

  testthat::expect_true(nrow(localFileRes_both)>0)

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAS_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "DELTAS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  # localFileRes <- read.table(fileToRead, sep="\t")

  testthat::expect_true(nrow(localFileRes_both)>0)


  # test deltaq creation
  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)

})

