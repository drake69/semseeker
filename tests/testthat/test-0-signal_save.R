test_that("signal-save",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  semseeker:::signal_save(signal_data,mySampleSheet,batch_id )

  signal_file <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE","")
  testthat::expect_true(file.exists(signal_file))


})
