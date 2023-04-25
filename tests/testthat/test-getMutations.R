testthat::test_that("mutations_get",{


  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- init_env(tempFolder)

  ####################################################################################

  mutations <- mutations_get(
               values = methylation_data[,1],
               figure = "HYPO",
               thresholds = tresholds,
               probe_features = probe_features,
               sampleName = mySampleSheet[1,"Sample_ID"]
               )

  expect_false(length(mutations)==0)

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  close_env()
})
