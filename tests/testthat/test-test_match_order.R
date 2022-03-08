test_that("test_match_order", {

  library(stringi)
  # tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  # envir <- init_env(tempFolder)
  #
  # nitem <- 5e4
  # nsamples <- 5
  # methylationData <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  # methylationData <- as.data.frame(matrix(methylationData,nitem,nsamples))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  # row.names(methylationData) <- probeFeatures$PROBE

  probeFeatures$ABSOLUTE <- paste(probeFeatures$CHR, probeFeatures$START, sep="_")

  #same order
  expect_true( test_match_order( probeFeatures$ABSOLUTE,probeFeatures$ABSOLUTE  ) )

  #ordre not matching
  expect_true( !test_match_order( probeFeatures$ABSOLUTE,sort(probeFeatures$ABSOLUTE, decreasing = TRUE)))

  #values not matching
  expect_true( !test_match_order( probeFeatures[-nrow(probeFeatures),"ABSOLUTE"],probeFeatures[-1, "ABSOLUTE"] ))

}
)
