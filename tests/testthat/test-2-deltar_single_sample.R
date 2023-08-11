test_that("deltar_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  # for (s in 1:ncol(methylation_data))
  {
    s <- 2
    message(s)
    sample_id <- colnames(methylation_data)[s]
    sample_detail <- mySampleSheet[ mySampleSheet$Sample_ID==sample_id ,c("Sample_ID","Sample_Group")]

    semseeker:::deltar_single_sample(
      values = methylation_data[,sample_id],
      high_thresholds = beta_superior_thresholds,
      low_thresholds = beta_inferior_thresholds,
      sample_detail = sample_detail,
      probe_features = probe_features
    )

    semseeker:::delta_single_sample(
      values = methylation_data[,sample_id],
      high_thresholds = beta_superior_thresholds,
      low_thresholds = beta_inferior_thresholds,
      sample_detail = sample_detail,
      probe_features = probe_features
    )

    semseeker:::analyze_single_sample( values = methylation_data[,sample_id],
      sliding_window_size = sliding_window_size,
      thresholds = beta_superior_thresholds ,
      figure = "HYPER",
      sample_detail = sample_detail ,
      bonferroni_threshold = 0.05,
      probe_features = probe_features)

    semseeker:::analyze_single_sample( values = methylation_data[,sample_id],
      sliding_window_size = sliding_window_size,
      thresholds = beta_inferior_thresholds ,
      figure = "HYPO",
      sample_detail = sample_detail ,
      bonferroni_threshold = 0.05,
      probe_features = probe_features)

    semseeker:::analyze_single_sample_both(sample_detail,"MUTATIONS")

    result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"MUTATIONS_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"MUTATIONS","BOTH"), "bed")
    mutations <- read.table(fileName)


    result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAR_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAR","BOTH"), "bedgraph")
    deltar <- read.table(fileName)
    # deltar <- read.csv(fileName, sep="\t")
    testthat::expect_true(nrow(deltar)>0)
    testthat::expect_true(nrow(deltar)==nrow(mutations))

    ####################################################################################

    outputFolder <- semseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAS_BOTH"))
    fileName <- semseeker:::file_path_build(outputFolder,c(sample_id,"DELTAS","BOTH"), "bedgraph")
    deltas <- read.table(fileName)
    testthat::expect_true(nrow(deltas)>0)
    testthat::expect_true(nrow(deltar)==nrow(deltas))

    ####################################################################################

    testthat::expect_true(nrow(mutations)==nrow(deltas))
  }

  ####################################################################################
  unlink(tempFolder, recursive = TRUE)
  semseeker:::close_env()

})
