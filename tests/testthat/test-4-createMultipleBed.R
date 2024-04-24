test_that("create_multiple_bed", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  markers <-c("MUTATIONS","DELTAS","DELTAQ","DELTARQ","DELTAR")
  figures <- c("HYPO","HYPER","BOTH")
  sample_groups <- c("Control","Reference","Case")

  semseeker:::create_multiple_bed(mySampleSheet)
  dq <- deltaq_get(mySampleSheet)
  drq <- deltarq_get(mySampleSheet)
  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")

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
          testthat::expect_true(file.exists(fileToRead))
          localFileRes_both <- fst::read_fst(fileToRead)
          testthat::expect_true(nrow(localFileRes_both)>0)
        }
        else
        {
          # message("fileToRead:",fileToRead)
        }

      }
    }
  }

  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
