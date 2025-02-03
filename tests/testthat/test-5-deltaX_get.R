test_that("deltaX_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,bonferroni_threshold = bonferroni_threshold)

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.01


  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  semseeker:::create_multiple_bed(mySampleSheet)
  resultPopulation <- semseeker:::deltaq_get(mySampleSheet)
  resultPopulation <- semseeker:::deltarq_get(mySampleSheet)
  resultPopulation <- semseeker:::deltap_get(mySampleSheet)
  resultPopulation <- semseeker:::deltarp_get(mySampleSheet)

  markers <-c("DELTAP","DELTAQ","DELTARP","DELTARQ")
  figures <- c("HYPO","HYPER","BOTH")
  sample_groups <- c("Control","Reference","Case")

  result_folderData  <-  semseeker:::dir_check_and_create(tempFolder, "Data")

  mutations_bed <- data.frame()
  marker_multiple_bed <- data.frame()
  count_multiple_bed <- 0
  total_sum <- data.frame()
  for (sample_group in sample_groups)
  {
    # sample_group <- "Case"
    for (marker in markers)
    {
      # marker <- "DELTARP"
      for (figure in figures)
      {
        # figure <- "HYPO"
        tempresult_folder <- semseeker:::dir_check_and_create(result_folderData,c(sample_group,paste(marker,figure, sep="_")))
        fileToRead <- semseeker:::file_path_build(tempresult_folder, c("MULTIPLE", marker ,figure ), "fst")
        if (file.exists(fileToRead))
        {
          mutations_bed <- plyr::rbind.fill(mutations_bed, semseeker:::read_multiple_bed (figure = figure, marker = "MUTATIONS", sample_group = sample_group))
          marker_multiple_bed <- plyr::rbind.fill(marker_multiple_bed, semseeker:::read_multiple_bed(figure = figure, marker = marker, sample_group = sample_group))
          # message(nrow(mutations_bed)==nrow(marker_multiple_bed))
          count_multiple_bed <- count_multiple_bed + 1
        }
      }
      total_sum <- plyr::rbind.fill(total_sum, data.frame("MARKER"=marker,"SUM"=sum(marker_multiple_bed[,4])))
      testthat::expect_true(nrow(marker_multiple_bed)==nrow(mutations_bed))
      testthat::expect_true(max(marker_multiple_bed[,4])==as.numeric(ssEnv$DELTAP_B))
      testthat::expect_true(min(marker_multiple_bed[,4])==1)
    }
  }

  # verify the markers produce a total burden different from each other
  testthat::expect_true(nrow(total_sum)==nrow(unique(total_sum$SUM)))
  testthat::expect_true(count_multiple_bed>0)

  semseeker:::close_env()
})

