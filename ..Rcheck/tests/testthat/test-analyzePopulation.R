  test_that("analize_population", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e3
  nsamples <- 20

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(beta_superior_thresholds) <- probe_features$PROBE
  row.names(beta_inferior_thresholds) <- probe_features$PROBE
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sample_sheet <- data.frame(Sample_Group, Sample_ID)
  beta_medians <- beta_superior_thresholds + beta_inferior_thresholds / 2
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  # browser()
  sp <- analize_population(envir = envir,
    methylation_data=methylation_data,
                    sliding_window_size = sliding_window_size,
                    beta_superior_thresholds = beta_superior_thresholds,
                    beta_inferior_thresholds = beta_inferior_thresholds,
                    sample_sheet = sample_sheet,
                    beta_medians = beta_medians,
                    bonferroni_threshold = bonferroni_threshold,
                    probe_features = probe_features
                    )

  sp$Sample_Group <- sample_sheet$Sample_Group

  message(nrow(sp))
  message(nrow(sample_sheet))
  expect_true(nrow(sp)==nrow(sample_sheet))

  future::plan( future::multisession)

})

