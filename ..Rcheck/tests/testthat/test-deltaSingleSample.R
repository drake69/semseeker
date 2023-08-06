testthat::test_that("delta_single_sample",{

  ssEnv <- init_env(tempFolder)

  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")

  nitem <- 1e4
  values <- as.data.frame(rnorm(nitem, mean=0.5, sd=0.7))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  high_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  low_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(values) <- probe_features$PROBE
  row.names(high_thresholds) <- probe_features$PROBE
  row.names(low_thresholds) <- probe_features$PROBE

  beta_medians <- high_thresholds - low_thresholds

  delta_single_sample(

    values = values,
    high_thresholds = high_thresholds,
    low_thresholds = low_thresholds,
    sample_detail = data.frame("Sample_ID"= Sample_ID, "Sample_Group"="Control"),
    beta_medians = beta_medians,
    probe_features = probe_features
  )

  result_folderData  <-  dir_check_and_create(tempFolder, "Data")
  outputFolder <- dir_check_and_create(result_folderData,c("Control","DELTAS_BOTH"))
  fileName <- file_path_build(outputFolder,c(Sample_ID,"DELTAS","BOTH"), "bedgraph")
  testthat::expect_true(file.exists(fileName))


})
