test_that("annotate_bed", {

  # unlink(tempFolder, recursive = T)
  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, showprogress= TRUE)

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

  # semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)
  semseeker:::annotate_bed ()
  dir_bed <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c("Annotated"))

  marker <- "MUTATIONS"
  figure <- "HYPO"
  area <- c("GENE")
  subarea <- c("BODY")
  final_bed <- semseeker:::read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true( nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))
  # testthat::expect_true( nrow(final_bed$SAMPLE_GROUP)>1)

  area <- c("GENE")
  subarea <- c("WHOLE")
  final_bed <- read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)>0)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  area =  "GENE"
  subarea <- "BODY"
  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90)

  bedFileName <- semseeker:::file_path_build(dir_bed , c(area, "Annotated"),"fst")
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  area =  "GENE"
  subarea = "TSS1500"
  final_bed <- read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)>0)

  area =  "GENE"
  subarea = "WHOLE"
  markers <- c("DELTAQ")
  final_bed <- semseeker:::read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)>0)

  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  final_bed <- semseeker:::read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)==nrow(unique(final_bed)))

  final_bed <- semseeker:::read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)>0)

  markers <- c("DELTAR")
  final_bed <- semseeker:::read_annotated_bed(figure,marker,area,subarea)
  testthat::expect_true(nrow(final_bed)>0)

  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
