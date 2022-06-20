test_that("annotate_bed", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  figures <- c( "BOTH")
  anomalies <- c("DELTAS")
  metaareas <- c("GENE")

  envir <- init_env(result_folder =  tempFolder, parallel_strategy = "sequential", maxResources = 90, figures = "BOTH", anomalies = "DELTAS", metaareas = "GENE")

  nitem <- 5e4
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

  sp <- analize_population(envir,
                          methylation_data=methylation_data,
                          sliding_window_size = 11,
                          beta_superior_thresholds = beta_superior_thresholds,
                          beta_inferior_thresholds = beta_inferior_thresholds,
                          sample_sheet = sample_sheet,
                          beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                          bonferroni_threshold = 0.01,
                          probe_features = probe_features
  )

  create_multiple_bed(envir, sample_sheet = sample_sheet, sp)

  populations <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS","DELTAS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  bedFileName <- file_path_build(envir$result_folderData , c(columnLabel, "ANNOTATED"),"csv")

  # file extsits
  expect_true(file.exists(bedFileName))

  # not empty data set
  expect_true(nrow(final_bed)>0)

  # has the correct header
  expect_true( columnLabel %in% colnames(final_bed))

  #read again  existent
  final_bed <- annotate_bed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  expect_true( columnLabel %in% colnames(final_bed))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)


  groups <- c("")
  probes_prefix = "PROBES"
  columnLabel =  "CHR"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  bedFileName <- file_path_build(envir$result_folderData , c(columnLabel, "ANNOTATED"),"csv")

})
