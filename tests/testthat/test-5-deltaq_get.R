test_that("deltaq_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  # browser()
  sp <- semseeker:::analyze_population(methylation_data=methylation_data,

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

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_HYPER")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"HYPER" ), "fst")
  deltaq_localFileRes_hyper <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(deltaq_localFileRes_hyper$VALUE))==0)
  testthat::expect_true(nrow(deltaq_localFileRes_hyper)>0)
  ####################################################################################

  tempresult_folder <- file.path(tempFolder,"Data","Control","MUTATIONS_HYPER")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"HYPER" ), "fst")
  mut_localFileRes_hyper <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(mut_localFileRes_hyper)==nrow(deltaq_localFileRes_hyper))

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAS_HYPER")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAS" ,"HYPER" ), "fst")
  deltas_localFileRes_hyper <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(mut_localFileRes_hyper)==nrow(deltas_localFileRes_hyper))

  testthat::expect_true(nrow(deltaq_localFileRes_hyper)==nrow(deltas_localFileRes_hyper))

  semseeker:::close_env()
})

