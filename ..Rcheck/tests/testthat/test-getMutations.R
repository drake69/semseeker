testthat::test_that("mutations_get",{

  init_env(tempFolder)

  Sample_ID <- stringi::stri_rand_strings(1, 15, pattern = "[A-Za-z]")

  nitem <- 1e4
  thresholds <- data.frame("thresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  row.names(thresholds) <- probe_features$PROBE
  row.names(values) <- row.names(thresholds)

  mutations <- mutations_get(
               values = values,
               figure = "HYPO",
               thresholds = thresholds,
               probe_features = probe_features,
               sampleName = Sample_ID
               )

  expect_false(length(mutations)==0)
})
