test_that("sortByCHRandSTART", {

  nitem <- 5e2

  probe_features <- PROBES[!is.na(PROBES$CHR),]
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- sortByCHRandSTART( probe_features[order(probe_features$START),])

  expect_true( test_match_order( sortByCHRandSTART(probe_features)$ABSOLUTE,second$ABSOLUTE))
}
)
