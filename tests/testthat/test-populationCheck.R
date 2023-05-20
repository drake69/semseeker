test_that(" semseeker:::population_check", {


  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  semseeker:::init_env(tempFolder)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  # all fine
  testthat::expect_true(is.null( semseeker:::population_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  testthat::expect_true(!is.null( semseeker:::population_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  testthat::expect_true(!is.null( semseeker:::population_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  testthat::expect_true(!is.null( semseeker:::population_check(mySampleSheet, methylation_data)))

  ####################################################################################

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

