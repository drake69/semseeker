test_that("coverage-analysis", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  semseeker:::get_meth_tech(methylation_data)
  res_cov <- semseeker:::coverage_analysis(methylation_data)
  testthat::expect_true(nrow(res_cov)!=0)

  # test cases of some ara fully missed
  nitem <- 1e2
  nsamples <- 21

  probe_features <- semseeker::pp_tot
  probe_features <- probe_features[!is.na(probe_features$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  row.names(methylation_data) <- probe_features$PROBE
  ####################################################################################

  res_cov <- semseeker:::coverage_analysis(methylation_data)
  testthat::expect_true(nrow(res_cov)!=0)

  ####################################################################################

  # semseeker:::close_env()
})

