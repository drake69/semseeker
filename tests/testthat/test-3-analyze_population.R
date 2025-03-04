test_that("analyze_population", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  mySampleSheet <- mySampleSheet[mySampleSheet$Sample_Group!="Reference",]
  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  testthat::expect_true(nrow(sp)==nrow(mySampleSheet))
  testthat::expect_true(all(mySampleSheet$Sample_ID %in% sp$Sample_ID))

  # check the column DELTAS_HYPO exists
  testthat::expect_true("DELTAS_HYPO" %in% colnames(sp))
  testthat::expect_true("DELTAS_HYPER" %in% colnames(sp))

  testthat::expect_true("DELTAR_HYPO" %in% colnames(sp))
  testthat::expect_true("DELTAR_HYPER" %in% colnames(sp))

  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_HYPER"])>0)>0)
  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_HYPO"])>0)>0)


  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})

