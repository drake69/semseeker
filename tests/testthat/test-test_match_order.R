test_that("test_match_order", {

  ####################################################################################
  #same order
  testthat::expect_true( test_match_order( probe_features$ABSOLUTE,probe_features$ABSOLUTE  ) )

  ####################################################################################
  #ordre not matching
  testthat::expect_true( !test_match_order( probe_features$ABSOLUTE,sort(probe_features$ABSOLUTE, decreasing = TRUE)))

  ####################################################################################
  #values not matching
  testthat::expect_true( !test_match_order( probe_features[-nrow(probe_features),"ABSOLUTE"],probe_features[-1, "ABSOLUTE"] ))
  # close_env()
}
)
