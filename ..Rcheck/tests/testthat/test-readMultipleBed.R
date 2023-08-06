test_that("read_multiple_bed", {

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

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("MUTATIONS","LESIONS","DELTAS")

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"
  sample_group <- unique(Sample_Group)

  probe_features <- get(paste0(probes_prefix,"WHOLE",sep=""), envir=asNamespace("semseeker"))


  res <-read_multiple_bed ( "MUTATIONS", "BOTH", probe_features, area, sample_group, groupingColumnLabel)
  testthat::expect_true(nrow(res)>0)

  res <-read_multiple_bed ( "DELTAS", "BOTH", probe_features, area, sample_group, groupingColumnLabel)
  testthat::expect_true(nrow(res)>0)
  # res <-read_multiple_bed ( "DELTAS", "HYPO", probe_features, area, sample_group, groupingColumnLabel)
  # res <-read_multiple_bed ( "DELTAS", "HYPER", probe_features, area, sample_group, groupingColumnLabel)
  # testthat::expect_true(nrow(res)>0)

}
)
