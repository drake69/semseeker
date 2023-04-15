test_that("annotate_bed", {

  
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  figures <- c( "BOTH")
  anomalies <- c("DELTAS","DELTAQ")
  metaareas <- c("GENE")
  ssEnv <- init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "BOTH", anomalies = "DELTAS", metaareas = "GENE")


  ####################################################################################

  get_meth_tech(methylation_data)

  ####################################################################################
  sp <- analize_population(methylation_data=methylation_data,
                          sliding_window_size = 11,
                          beta_superior_thresholds = beta_superior_thresholds,
                          beta_inferior_thresholds = beta_inferior_thresholds,
                          sample_sheet = mySampleSheet,
                          beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                          bonferroni_threshold = 0.01,
                          probe_features = probe_features
  )
  # sp$Sample_Group <- mySampleSheet$Sample_Group

  create_multiple_bed( sample_sheet = mySampleSheet)

  populations <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("DELTAS")

  groups <- c("")
  probes_prefix = "PROBES"
  columnLabel =  "PROBE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  testthat::expect_true( nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  groups <- c("CHR")
  probes_prefix = "PROBES_CHR_"
  columnLabel =  "CHR"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  # testthat::expect_true( columnLabel %in% colnames(final_bed))
  testthat::expect_true( nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  bedFileName <- file_path_build(ssEnv$result_folderData , c(columnLabel, "ANNOTATED"),"fst")

  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  anomalies <- c("DELTAQ")
  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)


  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  # not empty data set
  testthat::expect_true(nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  # has the correct header
  testthat::expect_true( columnLabel %in% colnames(final_bed))

  #read again  existent
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  testthat::expect_true( columnLabel %in% colnames(final_bed))


  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  # bedFileName <- file_path_build(ssEnv$result_folderData , c(columnLabel, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)

  # testthat::expect_true( columnLabel %in% colnames(final_bed))
  # bedFileName <- file_path_build(ssEnv$result_folderData , c(columnLabel, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)
  unlink(tempFolder,recursive = TRUE)
  close_env()
})
