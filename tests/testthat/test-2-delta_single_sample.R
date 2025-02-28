test_that("delta_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, inpute="median")

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  dss <- delta_single_sample(
    values = signal_data[,1],
    high_thresholds = signal_thresholds$signal_superior_thresholds,
    low_thresholds = signal_thresholds$signal_inferior_thresholds,
    sample_detail = mySampleSheet[1,c("Sample_ID","Sample_Group")],
    signal_medians = signal_medians,
    probe_features = probe_features
  )

  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
  outputFolder <- semseeker:::dir_check_and_create(result_folderData,c("Control","DELTAS_HYPER"))
  fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAS","HYPER"), "bedgraph", add_gz = TRUE)
  testthat::expect_true(file.exists(fileName))

  # message("fileName: ", fileName)
  # test I can open it
  res <- read.table(gzfile(fileName), header = FALSE)
  # message("res: ", res)
  testthat::expect_true(nrow(res)>0)

  ####################################################################################
  # result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
  # outputFolder <- semseeker:::dir_check_and_create(result_folderData,c("Control","DELTAS_HYPERS"))
  # fileName <- semseeker:::file_path_build(outputFolder,c(mySampleSheet[1,c("Sample_ID")],"DELTAS","HYPERS"), "bedgraph", add_gz = TRUE)
  # testthat::expect_true(file.exists(fileName))
  #
  # # message("fileName: ", fileName)
  # # test I can open it
  # res <- read.table(gzfile(fileName), header = FALSE)
  # # message("res: ", res)
  # testthat::expect_true(nrow(res)> 0)

  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)

})
