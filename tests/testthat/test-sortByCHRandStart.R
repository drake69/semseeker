test_that("sort_by_chr_and_start", {

  nitem <- 1e3

  probes <- probes_get("PROBES_Gene_","Whole")
  probe_features <- probes[!is.na(probes$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- sort_by_chr_and_start( probe_features[order(probe_features$START),])

  expect_true( test_match_order( sort_by_chr_and_start(probe_features)$ABSOLUTE,second$ABSOLUTE))

  ####################################################################################

  # close_env(envir)
}
)
