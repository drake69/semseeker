testthat::test_that("mutations_get",{

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 15, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder)

  Sample_ID <- stringi::stri_rand_strings(1, 15, pattern = "[A-Za-z]")

  nitem <- 4e5
  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
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

  expect_false(length(mutations)==0)
})
