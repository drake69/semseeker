test_that("test_match_order", {

  library(stringi)

  nitem <- 5e2
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
