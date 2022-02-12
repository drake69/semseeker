testthat::test_that("getMutations",{

  library(stringi)
  # tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 15, pattern = "[A-Za-z0-9]"),sep="")
  # init_env(tempFolder)

  Sample_ID <- stri_rand_strings(1, 15, pattern = "[A-Za-z]")

  nitem <- 5e4
  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  row.names(tresholds) <- probeFeatures$PROBE
  row.names(values) <- row.names(tresholds)

  mutations <- getMutations(
               values = values,
               figure = "HYPO",
               thresholds = tresholds,
               probeFeatures = probeFeatures,
               sampleName = Sample_ID
               )

  expect_false(length(mutations)==0)
})
