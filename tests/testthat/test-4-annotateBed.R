test_that("annotate_bed", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  # figures <- c("BOTH")
  # markers <- c("DELTAS","DELTAQ","DELTAR","DELTARQ")
  # markers <- c("BETA")
  # area <- c("GENE","CHR")
  # subarea <- c("BODY")
  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sp <- semseeker:::analyze_batch(methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    iqrTimes = 3,
    batch_id = 1
  )
  # sp$Sample_Group <- mySampleSheet$Sample_Group

  semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)

  sample_groups <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("DELTAS")
  area <-  c("CHR")
  subarea <- "WHOLE"
  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea
  )

  # testthat::expect_true( area %in% colnames(final_bed))
  testthat::expect_true( nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

    area <- c("GENE")
  subarea <- c("BODY")

  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures,
    markers,
    area,
    subarea
  )

  testthat::expect_true(nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))


  # subareas <- c("BODY")
  area =  "GENE"
  # ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, subareas = subareas)

  subarea <- "BODY"
  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)

  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  old_nrow <- nrow(unique(final_bed))
  # testing when added a subarea
  # subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  # ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, subareas = subareas)
  area =  "GENE"
  subarea = "TSS1500"
  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)

  bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area,subarea, "ANNOTATED"),"fst")
  testthat::expect_true(nrow(final_bed)!=old_nrow)

  area =  "GENE"
  subarea = "WHOLE"
  markers <- c("DELTAQ")
  # ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, markers = markers)
  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)

  # file extsits
  testthat::expect_true(nrow(final_bed)>0)

  # not empty data set
  # testthat::expect_true(nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  # has the correct header
  # testthat::expect_true( area %in% colnames(final_bed))

  #read again  existent
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)

  # testthat::expect_true( area %in% colnames(final_bed))
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  markers <- c("DELTARQ")
  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)


  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  markers <- c("DELTAR")
  # ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, markers = markers)
  # create and read
  bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area,subarea, "ANNOTATED"),"fst")
  final_bed <- semseeker:::annotate_bed (
    sample_groups ,
    figures ,
    markers ,
    area,
    subarea)

  # file extsits
  testthat::expect_true(file.exists(bedFileName))

  bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area,subarea, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)

  # testthat::expect_true( area %in% colnames(final_bed))
  # bedFileName <- semseeker:::file_path_build(ssEnv$result_folderData , c(area, "ANNOTATED"),"fst")
  # tt <- fst::read.fst(bedFileName)
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
