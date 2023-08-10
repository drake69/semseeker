test_that("create_heatmap", {

  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################


  batch_id <- 1
  sp <- semseeker:::analyze_batch(methylation_data =  methylation_data,
    sample_sheet =  mySampleSheet,
    sliding_window_size = sliding_window_size,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    batch_id = batch_id
  )

  semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)
  sp$Sample_Group <- mySampleSheet$Sample_Group


  semseeker:::create_excel_pivot()
  semseeker:::create_heatmap()

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/CHR/Control_Vs_Case_CHR_MUTATIONS_BOTH.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_MUTATIONS_BOTH.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_DELTAS_BOTH.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE/Control_Vs_Case_GENE_DELTAS_BOTH.png")))

  # unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
