test_that("test_match_order", {

  library(stringi)

  nitem <- 1e4

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
  # row.names(methylation_data) <- probe_features$PROBE

  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #same order
  testthat::expect_true( test_match_order( probe_features$ABSOLUTE,probe_features$ABSOLUTE  ) )

  #ordre not matching
  testthat::expect_true( !test_match_order( probe_features$ABSOLUTE,sort(probe_features$ABSOLUTE, decreasing = TRUE)))

  #values not matching
  testthat::expect_true( !test_match_order( probe_features[-nrow(probe_features),"ABSOLUTE"],probe_features[-1, "ABSOLUTE"] ))

}
)
