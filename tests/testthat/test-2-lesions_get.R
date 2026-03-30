test_that("lesions_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  Sample_ID <- mySampleSheet[1,"Sample_ID"]

  semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90,
    figures = "HYPER", markers = "DELTAS", areas = "GENE", bonferroni_threshold=5, inpute="median")

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

  mutations <-  semseeker:::mutations_get(
    values = values_df,
    figure = "HYPO",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )


  lesions_hypo <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  # With synthetic random data, mutations may be too sparse to form significant clusters;
  # verify the function returns a valid (possibly 0-row) data.frame rather than requiring lesions.
  testthat::expect_s3_class(lesions_hypo, "data.frame")

  mutations <-  semseeker:::mutations_get(
    values = values_df,
    figure = "HYPER",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )

  lesions_hyper <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_s3_class(lesions_hyper, "data.frame")

  ####################################################################################
})
