testthat::test_that("getLesions",{

  library(stringi)
  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")

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

  lesions <- getLesions(
    slidingWindowSize = 11,
    bonferroniThreshold = 0.05,
    mutationAnnotatedSorted = mutations,
    grouping_column = "CHR"
  )

  expect_false(length(lesions)==0)

})
