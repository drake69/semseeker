test_that("analyze_population", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]


  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  sp$Sample_Group <- mySampleSheet$Sample_Group

  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))


  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})

