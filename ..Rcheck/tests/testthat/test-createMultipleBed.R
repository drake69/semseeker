test_that("create_multiple_bed", {

  ssEnv <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e4
  nsamples <- 5
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(beta_superior_thresholds) <- probe_features$PROBE
  row.names(beta_inferior_thresholds) <- probe_features$PROBE
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sample_sheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analyze_population(methylation_data=methylation_data,
    sliding_window_size = 11,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features,
    bonferroni_threshold = 0.01,
  )

  sp$Sample_Group <- sample_sheet$Sample_Group
  create_multiple_bed( sample_sheet)

  tempresult_folder <-dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  # localFileRes_both <- read.table(fileToRead, sep="\t")
  localFileRes_both <- fst::read_fst(fileToRead)

  testthat::expect_true(nrow(localFileRes_both)>0)

  # tempresult_folder <-dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPO"))
  # fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"HYPO" ), "fst")
  # localFileRes_hypo <- read.table(fileToRead, sep="\t")
  #
  # tempresult_folder <-dir_check_and_create(ssEnv$result_folderData,c("Control","MUTATIONS_HYPER"))
  # fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"HYPER" ), "fst")
  # localFileRes_hyper <- read.table(fileToRead, sep="\t")

  # testthat::expect_true(nrow(localFileRes_hyper)>0 | nrow(localFileRes_hypo)>0 | nrow(localFileRes_both)>0)


})
