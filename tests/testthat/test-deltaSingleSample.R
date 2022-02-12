testthat::test_that("deltaSingleSample",{

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder)

  Sample_ID <- stri_rand_strings(1, 7, pattern = "[A-Za-z]")

  nitem <- 5e4
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.5, sd=0.7))

  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]

  highThresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  lowThresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(values) <- probeFeatures$PROBE
  row.names(highThresholds) <- probeFeatures$PROBE
  row.names(lowThresholds) <- probeFeatures$PROBE

  betaMedians <- highThresholds - lowThresholds

  deltaSingleSample(
    values = values,
    highThresholds = highThresholds,
    lowThresholds = lowThresholds,
    sampleDetail = data.frame("Sample_ID"= Sample_ID, "Sample_Group"="Control"),
    betaMedians = betaMedians,
    probeFeatures = probeFeatures
  )

  outputFolder <- dir_check_and_create(resultFolderData,c("Control","DELTAS_METHYLATION"))
  fileName <- file_path_build(outputFolder,c(Sample_ID,"DELTAS","METHYLATION"), "bedgraph")
  expect_true(file.exists(fileName))
})
