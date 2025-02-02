test_that("create_excel_pivot", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, start_fresh=TRUE,
    showprogress = FALSE, verbosity=verbosity)

  ####################################################################################
  # save sampèle sheet
  write.csv2(mySampleSheet, file.path(ssEnv$result_folderData,"sample_sheet_result.csv"), row.names = FALSE)

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  semseeker:::create_multiple_bed(mySampleSheet)
  dq <- deltaq_get(mySampleSheet)
  drq <- deltarq_get(mySampleSheet)
  dq <- deltap_get(mySampleSheet)
  drq <- deltarp_get(mySampleSheet)
  # semseeker:::annotate_bed()

  ####################################################################################

  # create and read
  semseeker:::create_excel_pivot()
  areas=c("GENE","ISLAND","DMR","CHR","PROBE", "SIGNAL")
  subarea=c("BODY","TSS1500","TSS200","1STEXON","3UTR","5UTR","EXONBND","WHOLE", "N_SHORE","S_SHORE","N_SHELF","S_SHELF")
  markers <-c("MUTATIONS","DELTAS","DELTAQ","DELTARQ","DELTAR")
  figures <- c("HYPO","HYPER","BOTH","MEAN","BOTHS")

  area = "GENE"
  subarea="WHOLE"
  marker="MUTATIONS"
  figure="BOTH"

  pivot_count <- 0
  for (area in areas)
  {
    for (subarea in subarea)
    {
      for (marker in markers)
      {
        for (figure in figures)
        {
          tempresult_folder <- semseeker:::dir_check_and_create(baseFolder = ssEnv$result_folderData , subFolders = c("Pivots",marker))
          fileToRead <- semseeker:::file_path_build(tempresult_folder, detailsFilename = c(marker,figure,area,subarea) , "csv", add_gz = TRUE)
          if(file.exists(fileToRead))
            pivot_count <- pivot_count + 1
        }
      }
    }
  }

  #at least  pivot
  testthat::expect_true(pivot_count > 0)

  semseeker:::close_env()
  unlink(tempFolder, recursive=TRUE)
})
