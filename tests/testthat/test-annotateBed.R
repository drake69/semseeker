test_that("annotate_bed", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  figures <- c( "BOTH")
  markers <- c("DELTAS","DELTAQ","DELTAR","DELTARQ")
  # markers <- c("BETA")
  areas <- c("GENE")
  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90)


  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sp <- semseeker:::analize_population(methylation_data=methylation_data,

    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features

  )
  # sp$Sample_Group <- mySampleSheet$Sample_Group

  semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)

  sample_groups <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("DELTAS")

  subareas <- c("")
  probes_prefix = "PROBES"
  area =  "PROBE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  testthat::expect_true(nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  subareas <- c("CHR")
  probes_prefix = "PROBES_CHR_"
  area =  "CHR"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  # testthat::expect_true( area %in% colnames(final_bed))
  testthat::expect_true( nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  subareas <- c("BODY")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  old_nrow <- nrow(unique(final_bed))
  # testing when added a subarea
  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  testthat::expect_true(nrow(final_bed)>old_nrow)

  markers <- c("DELTAQ")
  # create and read
  final_bed <- semseeker:::annotate_bed (

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
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  # has the correct header
  testthat::expect_true( area %in% colnames(final_bed))

  #read again  existent
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  testthat::expect_true( area %in% colnames(final_bed))
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  markers <- c("DELTARQ")
  # create and read
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)


  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  markers <- c("DELTAR")
  # create and read
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)


  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)

  # testthat::expect_true( area %in% colnames(final_bed))
  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
