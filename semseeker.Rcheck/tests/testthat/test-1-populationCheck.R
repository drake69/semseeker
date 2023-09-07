test_that(" semseeker:::sample_group_check", {


  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  # all fine
  testthat::expect_true(is.null( semseeker:::sample_group_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  testthat::expect_true(!is.null( semseeker:::sample_group_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  testthat::expect_true(!is.null( semseeker:::sample_group_check(mySampleSheet, methylation_data)))

  ####################################################################################

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  testthat::expect_true(!is.null( semseeker:::sample_group_check(mySampleSheet, methylation_data)))

  ####################################################################################

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})

