test_that("range_beta_values", {

  ssEnv <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e4
  nsamples <- 21

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  row.names(methylation_data) <- probe_features$PROBE

  # Sample_ID <- stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
  # colnames(methylation_data) <- Sample_ID
  # Sample_Group <- c(rep("Control",nsamples/3),rep("Case",nsamples/3),rep("Reference",nsamples/3))
  # sample_sheet <- data.frame(Sample_Group, Sample_ID)
  #
  # sliding_window_size <- 11
  # bonferroni_threshold <- 0.01
  # batch_id <- 1
  iqrTimes <- 3

  rr <- range_beta_values(methylation_data)
  # sp <- analyze_batch(
  #                       methylation_data =  methylation_data,
  #                       sample_sheet =  sample_sheet,
  #                       sliding_window_size = sliding_window_size,
  #                       bonferroni_threshold =  bonferroni_threshold,
  #                       iqrTimes =  iqrTimes,
  #                       batch_id = batch_id)
  #

  # message(nrow(sp))
  # message(nrow(sample_sheet))
  testthat::expect_true(sum(colnames(rr)==c("beta_inferior_thresholds","beta_superior_thresholds","beta_median_values"))==3)
  testthat::expect_true(nrow(rr)==nrow(methylation_data))
  testthat::expect_true(rr$beta_inferior_thresholds[1]!=rr$beta_superior_thresholds[1])
  testthat::expect_true(rr$beta_inferior_thresholds[1]!=rr$beta_median_values[1])
  # future::plan( future::multisession)
})

