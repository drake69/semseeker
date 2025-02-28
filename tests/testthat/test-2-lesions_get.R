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

  mutations <-  semseeker:::mutations_get(
    values = signal_data[,1],
    figure = "HYPO",
    thresholds = signal_thresholds$signal_inferior_thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )


  lesions_hypo <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hypo)!=0)

  mutations <-  semseeker:::mutations_get(
    values = signal_data[,1],
    figure = "HYPER",
    thresholds = signal_thresholds$signal_superior_thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions_hyper <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hyper)!=0)

  ####################################################################################
})
