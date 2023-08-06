test_that("sample_group_check", {


  ssEnv <- init_env(tempFolder)


  # all fine
  testthat::expect_true(is.null(sample_group_check(mySampleSheet, methylation_data)))

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  testthat::expect_true(!is.null(sample_group_check(mySampleSheet, methylation_data)))

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  testthat::expect_true(!is.null(sample_group_check(mySampleSheet, methylation_data)))

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  testthat::expect_true(!is.null(sample_group_check(mySampleSheet, methylation_data)))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)

})

