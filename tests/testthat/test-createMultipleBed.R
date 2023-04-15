test_that("create_multiple_bed", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  get_meth_tech(methylation_data)

  ####################################################################################

  sp <- analize_population(methylation_data=methylation_data,
                           sliding_window_size = 11,
                           beta_superior_thresholds = beta_superior_thresholds,
                           beta_inferior_thresholds = beta_inferior_thresholds,
                           sample_sheet = mySampleSheet,
                           beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                           bonferroni_threshold = 0.01,
                           probe_features = probe_features
  )

  create_multiple_bed(mySampleSheet)

  tempresult_folder <-dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(nrow(localFileRes_both)>0)

  ####################################################################################

  unlink(tempFolder, recursive = TRUE)
  close_env()
})
