testthat::test_that("lesions_get",{

  library(stringi)
  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")

  nitem <- 5e4
  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
  probe_features <- PROBES[!is.na(PROBES$CHR),]
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
  row.names(tresholds) <- probe_features$PROBE
  row.names(values) <- row.names(tresholds)

  mutations <- mutations_get(
    values = values,
    figure = "HYPO",
    thresholds = tresholds,
    probe_features = probe_features,
    sampleName = Sample_ID
  )

  lesions <- lesions_get(
    sliding_window_size = 11,
    bonferroni_threshold = 0.05,
    mutation_annotated_sorted = mutations,
    grouping_column = "CHR"
  )

  expect_false(length(lesions)==0)

})
