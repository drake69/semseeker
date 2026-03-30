test_that("analyze_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <-  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute = "mean")
  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  # Build a bed-like data.frame matching the production code path (values read from bed file)
  values_df <- data.frame(
    CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
    START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
    END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
    VALUE = as.numeric(signal_data[, 1])
  )

  sp <- semseeker:::analyze_single_sample(
    values = values_df,
    thresholds = signal_thresholds,
    figure = "HYPO",
    sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")]
  )


  outputFolder <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPO"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"MUTATIONS","HYPO"), "bed", add_gz = TRUE)
  testthat::expect_true(file.exists(fileName))

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})
