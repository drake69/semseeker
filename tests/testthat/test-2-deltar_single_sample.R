test_that("deltar_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(result_folder= tempFolder, bonferroni_threshold = 0.05, showprogress=showprogress, start_fresh = TRUE, inpute="median")

  ####################################################################################

  gg <- semseeker:::get_meth_tech(signal_data)
  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]
  sample_group <- "Reference"
  ####################################################################################
  # for (s in 1:ncol(signal_data))
  {
    s <- 2
    message(s)
    sample_id <- colnames(signal_data)[s]
    sample_detail <- mySampleSheet[ mySampleSheet$Sample_ID==sample_id & mySampleSheet$Sample_Group==sample_group ,c("Sample_ID","Sample_Group")]

    dr1 <- semseeker:::deltar_single_sample(
      values = signal_data[,sample_id],
      high_thresholds = signal_thresholds$signal_superior_thresholds,
      low_thresholds = signal_thresholds$signal_superior_thresholds,
      sample_detail = sample_detail,
      probe_features = probe_features
    )

    ds1 <- semseeker:::delta_single_sample(
      values = signal_data[,sample_id],
      high_thresholds = signal_thresholds$signal_superior_thresholds,
      low_thresholds = signal_thresholds$signal_superior_thresholds,
      sample_detail = sample_detail,
      probe_features = probe_features
    )

    ass <- semseeker:::analyze_single_sample( values = signal_data[,sample_id],
      thresholds = signal_thresholds$signal_superior_thresholds ,
      figure = "HYPER",
      sample_detail = sample_detail ,
      probe_features = probe_features)

    semseeker:::analyze_single_sample( values = signal_data[,sample_id],

      thresholds = signal_thresholds$signal_superior_thresholds ,
      figure = "HYPO",
      sample_detail = sample_detail ,
      probe_features = probe_features)

    mb <- semseeker:::analyze_single_sample_both(sample_detail,"MUTATIONS")

    result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"MUTATIONS_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"MUTATIONS","BOTH"), "bed", add_gz = TRUE)
    mutations <- read.table(fileName)


    result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAR_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAR","BOTH"), "bedgraph", add_gz = TRUE)
    deltar <- read.table(gzfile(fileName))
    # deltar <- read.csv(fileName, sep="\t")
    testthat::expect_true(nrow(deltar)>0)
    testthat::expect_true(nrow(deltar)==nrow(mutations))

    ####################################################################################

    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAS_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAS","BOTH"), "bedgraph", add_gz = TRUE)
    deltas <- read.table(fileName)
    testthat::expect_true(nrow(deltas)>0)
    testthat::expect_true(nrow(deltar)==nrow(deltas))

    ####################################################################################

    testthat::expect_true(nrow(mutations)==nrow(deltas))

    ####################################################################################

    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAR_HYPER"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAR","HYPER"), "bedgraph", add_gz = TRUE)
    deltar <- read.table(fileName)
    testthat::expect_true(nrow(deltar)>0)

    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAS_HYPER"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAS","HYPER"), "bedgraph", add_gz = TRUE)
    deltas <- read.table(fileName)

    testthat::expect_true(nrow(deltar)==nrow(deltas))

  }

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)

})
