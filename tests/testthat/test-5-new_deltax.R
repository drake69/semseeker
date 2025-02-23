test_that("deltaX_get", {

  tempFolder <- tempFolders[1]
  message(tempFolder)
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,bonferroni_threshold = bonferroni_threshold, inpute="median", start_fresh=FALSE)

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.05

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  keys <- ssEnv$keys_markers_figures_default
  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )
  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Case",],keys)

  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )
  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Control",],keys)

  for( k in 1:nrow(keys))
  {
    key <- keys[k,]
    file_name <- pivot_file_name(key)
  }

  mySampleSheet <- semseeker:::deltaX_get(mySampleSheet[mySampleSheet$Sample_Group!="Reference",])

  # check deltaX pivot are coherent with mutations
  # same row number
  # same column number
  # dataframe of deltaX>0 minus dataframe of mutations == 0

  semseeker:::annotate_position_pivots()

  semseeker:::close_env()
})

