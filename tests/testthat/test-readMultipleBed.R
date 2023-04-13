test_that("read_multiple_bed", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e3
  nsamples <- 5
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probes <- semseeker::PROBES
  probe_features <- probes[!is.na(probes$START),c("CHR","START","PROBE")]
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

  ####################################################################################

  get_meth_tech(methylation_data)

  ####################################################################################
  sp <- analize_population(methylation_data=methylation_data,
                           sliding_window_size = 11,
                           beta_superior_thresholds = beta_superior_thresholds,
                           beta_inferior_thresholds = beta_inferior_thresholds,
                           sample_sheet = sample_sheet,
                           beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                           bonferroni_threshold = 0.01,
                           probe_features = probe_features
  )

  sp$Sample_Group <- sample_sheet$Sample_Group
  create_multiple_bed(envir, sample_sheet)

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS","DELTAS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"
  populationName <- unique(Sample_Group)

  probe_features <- get(paste0(probes_prefix,"Whole",sep=""), envir = asNamespace("semseeker"))


  res <-read_multiple_bed (envir, "MUTATIONS", "BOTH", probe_features, columnLabel, populationName, groupingColumnLabel)
  expect_true(nrow(res)>0)

  res <-read_multiple_bed (envir, "DELTAS", "BOTH", probe_features, columnLabel, populationName, groupingColumnLabel)
  expect_true(nrow(res)>0)
  ####################################################################################

  # res <-read_multiple_bed (envir, "DELTAS", "HYPO", probe_features, columnLabel, populationName, groupingColumnLabel)
  # res <-read_multiple_bed (envir, "DELTAS", "HYPER", probe_features, columnLabel, populationName, groupingColumnLabel)
  # expect_true(nrow(res)>0)
  close_env(envir)
}
)
