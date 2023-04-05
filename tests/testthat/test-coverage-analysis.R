test_that("coverage-analysis", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e5
  nsamples <- 21

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  row.names(methylation_data) <- probe_features$PROBE

  res_cov <- coverage_analysis(methylation_data)

  # message(nrow(sp))
  # message(nrow(sample_sheet))
  testthat::expect_true(nrow(sp)==nrow(sample_sheet))
  testthat::expect_true(sum(na.omit(sp[,"MUTATIONS_BOTH"])>0)>0)

  close_env(envir)
})

