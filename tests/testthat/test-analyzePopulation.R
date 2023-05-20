test_that("analize_population", {



  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  semseeker:::init_env(tempFolder, parallel_strategy = "sequential")

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################


  # browser()
  sp <- semseeker:::analize_population(methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features
  )

  sp$Sample_Group <- mySampleSheet$Sample_Group

  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))

  ####################################################################################

  beta_file <- file_path_build(folder_to_save,c(mySampleSheet[1,"Sample_ID"],"BETAS",figure),"bed")
  testthat::expect_true(file.exists(beta_file))

  beta_file <- read.table(beta_file, sep="\t")
  testthat::expect_true(nrow(beta_file)==nrow(methylation_data))

  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()

})

