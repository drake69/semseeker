testthat::test_that("lesions_get",{


  Sample_ID <- mySampleSheet[1,"Sample_ID"]
  bonferroni_threshold <- 5

  mutations <-  semseeker:::mutations_get(
    values = methylation_data[,1],
    figure = "HYPO",
    thresholds = tresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions_hypo <- semseeker:::lesions_get(
    sliding_window_size = sliding_window_size,
    bonferroni_threshold = bonferroni_threshold,
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hypo)!=0)

  mutations <-  semseeker:::mutations_get(
    values = methylation_data[,1],
    figure = "HYPER",
    thresholds = tresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions_hyper <- semseeker:::lesions_get(
    sliding_window_size = sliding_window_size,
    bonferroni_threshold = bonferroni_threshold,
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hyper)!=0)

  ####################################################################################

})
