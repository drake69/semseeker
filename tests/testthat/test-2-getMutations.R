testthat::test_that(" semseeker:::mutations_get",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder)

  ####################################################################################

  mutations <-  semseeker:::mutations_get(
               values = methylation_data[,1],
               figure = "HYPO",
               thresholds = beta_inferior_thresholds,
               probe_features = probe_features,
               sampleName = mySampleSheet[1,"Sample_ID"]
               )

  expect_false(length(mutations)==0)

  ####################################################################################
  # expect the count of mutations is not less than 80% of total added mutations
  expect_true(sum(mutations$MUTATIONS==1)> 0.8 * perc_epimutation * nprobes/2 )

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()
})
