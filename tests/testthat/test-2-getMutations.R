test_that(" semseeker:::mutations_get",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(result_folder=tempFolder, inpute="median")

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  mutations <-  semseeker:::mutations_get(
               values = signal_data[,1],
               figure = "HYPO",
               thresholds = signal_thresholds$signal_inferior_thresholds,
               probe_features = probe_features,
               sampleName = mySampleSheet[1,"Sample_ID"]
               )

  expect_false(length(mutations)==0)

  ####################################################################################
  # expect the count of mutations is not less than 80% of total added mutations
  # expect_true(sum(mutations$MUTATIONS==1)> 0 )

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
