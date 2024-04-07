test_that("semeeker", {

  ssEnv <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e4
  nsamples <- 20

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  signal_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  signal_data <- as.data.frame(matrix(signal_data,nitem,nsamples))

  signal_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  signal_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(signal_superior_thresholds) <- probe_features$PROBE
  row.names(signal_inferior_thresholds) <- probe_features$PROBE
  row.names(signal_data) <- probe_features$PROBE

  Sample_ID <- stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
  colnames(signal_data) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sample_sheet <- data.frame(Sample_Group, Sample_ID)
  signal_medians <- signal_superior_thresholds + signal_inferior_thresholds / 2
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  # browser()
  sp <- analyze_population(signal_data=signal_data,
    sliding_window_size = 11,
    sliding_window_size = sliding_window_size,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features,
    bonferroni_threshold = 0.01,
  )

  create_multiple_bed( sample_sheet)

  resiltPopulation <- deltaq_get( sp)

  # test deltaq creation
  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)

})

