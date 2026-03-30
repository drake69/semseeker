test_that("lesions_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  Sample_ID <- mySampleSheet[1, "Sample_ID"]

  semseeker:::init_env(
    result_folder       = tempFolder,
    parallel_strategy   = parallel_strategy,
    maxResources        = 90,
    figures             = "HYPER",
    markers             = "DELTAS",
    areas               = "GENE",
    bonferroni_threshold = 5,
    inpute              = "median"
  )

  if (!exists("signal_thresholds")) {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data), ]

  values_df <- data.frame(
    CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
    START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
    END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
    VALUE = as.numeric(signal_data[, 1])
  )

  mutations <- semseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )

  # ── HYPO lesions: valid data.frame returned ────────────────────────────────
  lesions_hypo <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column           = "CHR"
  )

  testthat::expect_s3_class(lesions_hypo, "data.frame")

  # output columns must be exactly CHR / START / END
  testthat::expect_true(all(c("CHR", "START", "END") %in% colnames(lesions_hypo)))

  # lesion count cannot exceed probe count
  testthat::expect_true(nrow(lesions_hypo) <= nrow(mutations))

  # ── HYPER lesions ──────────────────────────────────────────────────────────
  mutations_hyper <- semseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )

  lesions_hyper <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations_hyper,
    grouping_column           = "CHR"
  )

  testthat::expect_s3_class(lesions_hyper, "data.frame")

  # ── Edge case: NULL input returns NULL ─────────────────────────────────────
  lesions_null <- semseeker:::lesions_get(
    mutation_annotated_sorted = NULL,
    grouping_column           = "CHR"
  )
  testthat::expect_null(lesions_null)

  # ── Edge case: empty data.frame returns 0-row result ──────────────────────
  lesions_empty <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations[0, ],
    grouping_column           = "CHR"
  )
  testthat::expect_equal(nrow(lesions_empty), 0)

  # ── Alternative grouping: PROBE column ────────────────────────────────────
  mutations_with_probe <- mutations
  mutations_with_probe$PROBE <- paste0("cg", formatC(seq_len(nrow(mutations)), width = 7, flag = "0"))
  lesions_probe <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations_with_probe,
    grouping_column           = "PROBE"
  )
  testthat::expect_s3_class(lesions_probe, "data.frame")

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
