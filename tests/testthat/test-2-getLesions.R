test_that("lesions_get",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  Sample_ID <- mySampleSheet[1,"Sample_ID"]

  semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "BOTH", markers = "DELTAS", areas = "GENE", bonferroni_threshold=5)

  mutations <-  semseeker:::mutations_get(
    values = signal_data[,1],
    figure = "HYPO",
    thresholds = signal_inferior_thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  # set to 1 10 values of mutations
  mutations[1:100,"MUTATIONS"] <- 1

  lesions_hypo <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hypo)!=0)

  mutations <-  semseeker:::mutations_get(
    values = signal_data[,1],
    figure = "HYPER",
    thresholds = signal_superior_thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  mutations[1:100,"MUTATIONS"] <- 1
  lesions_hyper <- semseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hyper)!=0)

  ####################################################################################
})
