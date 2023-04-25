test_that("sort_by_chr_and_start", {


  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- sort_by_chr_and_start( probe_features[order(probe_features$START),])

  testthat::expect_true( test_match_order( sort_by_chr_and_start(probe_features)$ABSOLUTE,second$ABSOLUTE))

  ####################################################################################

  # close_env()
}
)
