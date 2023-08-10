test_that("annotate_bed", {


  figures <- c( "BOTH")
  markers <- c("DELTAS","DELTAQ")
  areas <- c("GENE")

  init_env(result_folder =  tempFolder, parallel_strategy = "sequential", maxResources = 90, figures = "BOTH", markers = "DELTAS", areas = "GENE")

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

  create_multiple_bed( sample_sheet = sample_sheet)

  sample_groups <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("DELTAS")

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  bedFileName <- file_path_build(ssEnv$result_folderData , c(area, "Annotated"),"fst")


  markers <- c("DELTAQ")
  # create and read
  final_bed <- annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)


  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  # not empty data set
  testthat::expect_true(nrow(final_bed)>0)

  # has the correct header
  testthat::expect_true( area %in% colnames(final_bed))

  #read again  existent
  final_bed <- annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  testthat::expect_true( area %in% colnames(final_bed))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)


  subareas <- c("CHR")
  probes_prefix = "PROBES_CHR_"
  area =  "CHR"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  # testthat::expect_true( area %in% colnames(final_bed))
  testthat::expect_true( nrow(final_bed)>0)

  # bedFileName <- file_path_build(ssEnv$result_folderData , c(area, "Annotated"),"fst")
  # tt <- fst::read.fst(bedFileName)

  subareas <- c("")
  probes_prefix = "PROBES"
  area =  "PROBE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  testthat::expect_true( nrow(final_bed)>0)
  # testthat::expect_true( area %in% colnames(final_bed))
  # bedFileName <- file_path_build(ssEnv$result_folderData , c(area, "Annotated"),"fst")
  # tt <- fst::read.fst(bedFileName)

})
