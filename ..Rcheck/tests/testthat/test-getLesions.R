testthat::test_that("lesions_get",{

  library(stringi)
  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")

  nitem <- 1e4

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]


  thresholds <- data.frame("thresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))

  row.names(thresholds) <- probe_features$PROBE
  row.names(values) <- row.names(thresholds)

  mutations <- mutations_get(
    values = values,
    figure = "HYPO",
    thresholds = thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions_hypo <- lesions_get(
    sliding_window_size = 11,
    bonferroni_threshold = 5,
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  mutations <- mutations_get(
    values = values,
    figure = "HYPER",
    thresholds = thresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions_hyper <- lesions_get(
    sliding_window_size = 11,
    bonferroni_threshold = 5,
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  testthat::expect_true(nrow(lesions_hyper)!=0 | nrow(lesions_hypo)!=0)

})
