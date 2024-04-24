test_that("association_analysis", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  ####################################################################################

  semseeker:::semseeker( sample_sheet =  mySampleSheet,signal_data =  signal_data,result_folder = tempFolder,parallel_strategy=parallel_strategy, figures="BOTH",
    markers=c("DELTAQ"), areas=c("PROBE","CHR","GENE"), start_fresh = TRUE)

  ####################################################################################
  ssEnv <- get_session_info(tempFolder)
  markers <- c("DELTAQ")
  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",markers=markers, areas="GENE", subareas=c("WHOLE"))

  res <- get_results_inference(inference_details, markers[1], ssEnv$result_folderInference)
  inferenceFile <- semseeker:::inference_file_name(inference_details, markers[1], ssEnv$result_folderInference)
  testthat::expect_true(file.exists(inferenceFile))

  ####################################################################################


  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("poisson"),
    "transformation"="",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH", markers=markers, areas="PROBE")

  inferenceFile <- semseeker:::inference_file_name(inference_details, markers[1], ssEnv$result_folderInference)
  testthat::expect_true(file.exists(inferenceFile))

  ####################################################################################


  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg-permutation_0.5_5_10_0.9"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="PROBE")

  inferenceFile <- semseeker:::inference_file_name(inference_details, markers[1], ssEnv$result_folderInference)
  testthat::expect_true(file.exists(inferenceFile))

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Covariates1",
                                   "covariates"=c("Phenotest+Covariates2"),
                                   "family_test"=c("quantreg_0.5_5_10"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  # inference_details,result_folder, maxResources, parallel_strategy
  # test with area selection
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH",  markers=markers, areas="PROBE",areas_selection=c("chr1","chr2"))

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)
  firstRun <- nrow(localFileRes)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)
  secondRun <- nrow(localFileRes)

  ####################################################################################

  # inference_details,result_folder, maxResources, parallel_strategy
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH",
    markers=markers, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  thirdRun <- nrow(localFileRes)
  testthat::expect_true(thirdRun > firstRun)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  fourthRun <- nrow(localFileRes)
  testthat::expect_true(fourthRun > secondRun)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  # inference_details,result_folder, maxResources, parallel_strategy
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH", markers=markers,
    areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)


  ####################################################################################

  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("quantreg_0.5"),
  #   "transformation"="scale",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,
  #   figures="BOTH", markers=markers, areas="PROBE")
  #
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_5_10_0.99"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH", markers=markers, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("quantreg_0.5_5_10_0.99_np"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH", markers=markers, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  #todo: test incremental association analysis
  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_5_10_0.9"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="BOTH", markers=markers, areas="CHR")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="log",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################


  inference_details <- expand.grid("independent_variable"= "Group",
    "covariates"="",
    "family_test"=c("t.test"),
    "transformation"="log",
    "depth_analysis"=1,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Group",
    "covariates"="",
    "family_test"=c("pearson"),
    "transformation"="log",
    "depth_analysis"=1,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="quantile_5",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="quantile_5",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates" =c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="scale",
    "depth_analysis"=1,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="log10",
    "depth_analysis"=2,
    "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="logn",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)

  areas_selection <- rownames(signal_data)[1:100]
  # inference_details,result_folder, maxResources, parallel_strategy
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas_selection=areas_selection, areas="PROBE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################


  inference_details <- expand.grid("independent_variable"= "Sample_Group",
    "family_test"=c("mean-permutation_100_1000_0.95"),
    "transformation"="",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################


  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "family_test"=c("spearman-permutation_100_1000_0.95"),
    "transformation"="",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= c("Age","AgeRange5","AgeRange10","AgeRange20"),
      "family_test"=c("polynomial_4_1"),
      "covariates"="",
      "transformation"="scale",
      "depth_analysis"=2,
      "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

    ####################################################################################

  inference_details <- expand.grid("independent_variable"= c("Age","AgeRange5","AgeRange10","AgeRange20"),
      "family_test"=c("polynomial_4_1_predict"),
      "covariates"="",
      "transformation"="scale",
      "depth_analysis"=2,
      "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <- expand.grid("independent_variable"= c("Age","AgeRange5","AgeRange10","AgeRange20"),
      "family_test"=c("exp_1"),
      "covariates"="",
      "transformation"="scale",
      "depth_analysis"=2,
      "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################

  inference_details <-    expand.grid("independent_variable"= c("Age","AgeRange5","AgeRange10","AgeRange20"),
      "family_test"=c("log_1"),
      "covariates"="",
      "transformation"="scale",
      "depth_analysis"=2,
      "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="BOTH",
    markers=markers, areas="GENE")

  localFileRes <- get_results_inference(inference_details = inference_details,markers[1])
  testthat::expect_true(nrow(localFileRes)>0)

  localFileRes <- get_results_inference(inference_details = inference_details,markers[2])
  testthat::expect_true(nrow(localFileRes)>0)

  ####################################################################################


  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})
