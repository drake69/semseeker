test_that("deltarq_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  # browser()
  sp <- semseeker:::analyze_population(
    methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features
  )

  semseeker:::create_multiple_bed(mySampleSheet)
  resultPopulation <- deltarq_get(mySampleSheet)

  # test deltaq creation
  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTARQ_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTARQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)
  testthat::expect_true(nrow(localFileRes_both)>0)

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTARQ_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTARQ" ,"BOTH" ), "fst")
  deltarq_localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(deltarq_localFileRes_both$VALUE))==0)
  testthat::expect_true(nrow(deltarq_localFileRes_both)>0)
  ####################################################################################


  tempresult_folder <- file.path(tempFolder,"Data","Control","MUTATIONS_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  mut_localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(mut_localFileRes_both)>=nrow(deltarq_localFileRes_both))

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAS_BOTH")
  fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", "DELTAS" ,"BOTH" ), "fst")
  deltas_localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(deltas_localFileRes_both)>=nrow(deltarq_localFileRes_both))

  testthat::expect_true(nrow(mut_localFileRes_both)==nrow(deltas_localFileRes_both))
  ####################################################################################
  unlink(tempFolder, recursive = T)
  semseeker:::close_env()
})

