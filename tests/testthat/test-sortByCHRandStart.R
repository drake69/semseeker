test_that("sortByCHRandSTART", {

  nitem <- 5e2

  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]

  probeFeatures$ABSOLUTE <- paste(probeFeatures$CHR, probeFeatures$START, sep="_")

  #order not matching
  second <- sortByCHRandSTART( probeFeatures[order(probeFeatures$START),])

  expect_true( test_match_order( sortByCHRandSTART(probeFeatures)$ABSOLUTE,second$ABSOLUTE))
}
)
