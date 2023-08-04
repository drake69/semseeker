test_that("deltaq_get", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  # browser()
  sp <- semseeker:::analize_population(methylation_data=methylation_data,

    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features

  )

  semseeker:::create_multiple_bed(mySampleSheet)
  resultPopulation <- semseeker:::deltaq_get(mySampleSheet)

  # test deltaq creation
  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)
  testthat::expect_true(nrow(localFileRes_both)>0)

  tempresult_folder <- file.path(tempFolder,"Data","Reference","DELTAQ_HYPER")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"HYPER" ), "fst")
  localFileRes_hyper <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_hyper$VALUE))==0)
  testthat::expect_true(nrow(localFileRes_hyper)>0)
  ####################################################################################

  semseeker:::close_env()
})

