test_that("deltarq_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, bonferroni_threshold = bonferroni_threshold)

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01

  
  sp <- semseeker:::analyze_population(
    signal_data=signal_data,

    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  semseeker:::create_multiple_bed(mySampleSheet)
  resultPopulation <- deltarq_get(mySampleSheet)

  markers <-c("DELTARQ")
  figures <- c("HYPO","HYPER","BOTH")
  sample_groups <- c("Control","Reference","Case")

  semseeker:::create_multiple_bed(mySampleSheet)
  dq <- deltaq_get(mySampleSheet)
  drq <- deltarq_get(mySampleSheet)
  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")

  count_multiple_bed <- 0
  for (sample_group in sample_groups)
  {
    for (marker in markers)
    {
      for (figure in figures)
      {
        tempresult_folder <- semseeker:::dir_check_and_create(result_folderData,c(sample_group,paste(marker,figure, sep="_")))
        fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", marker ,figure ), "fst")
        if (file.exists(fileToRead))
        {
          read_multiple_bed <-semseeker:::read_multiple_bed (figure = figure, marker = marker, sample_group = sample_group)
          testthat::expect_true(nrow(read_multiple_bed)>0)
          count_multiple_bed <- count_multiple_bed + 1
        }
      }
    }
  }

  testthat::expect_true(count_multiple_bed>0)

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder, recursive = T)

})

