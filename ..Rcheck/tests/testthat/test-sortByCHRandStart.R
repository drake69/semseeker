test_that("sort_by_chr_and_start", {

  nitem <- 1e4

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- sort_by_chr_and_start( probe_features[order(probe_features$START),])

  testthat::expect_true( test_match_order( sort_by_chr_and_start(probe_features)$ABSOLUTE,second$ABSOLUTE))
}
)
