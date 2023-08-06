test_that("association_analysis", {


  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,parallel_strategy="sequential", figures="BOTH", markers="DELTAS", areas="PROBE")

  #todo: test incremental association analysis
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_1000_15000_0.9"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=1,
  #                                  "filter_p_value" = FALSE)
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential")
  #
  # fileToRead <- file_path_build(inferenceFolder, "1_Phenotest_scale_quantreg_0.5_1000_15000_0.9_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)


  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Phenotest+Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_1000_15000"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="PROBE")

  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("quantreg_0.5"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential",
    figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_100_1000_0.99"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("quantreg_0.5_100_1000_0.99_np"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="PROBE")


  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,parallel_strategy="sequential", figures="BOTH", markers="DELTAS", areas="CHR")
  inferenceFolder <- file.path(tempFolder,"Inference")


  #todo: test incremental association analysis

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_1000_15000_0.9"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", markers=c("DELTAS","DELTAQ"), areas="CHR")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_1000_15000_0.9_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  test_both <- nrow(localFileRes)
  testthat::expect_true(nrow(localFileRes)>0)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", markers="DELTAS", areas="CHR")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_1000_15000_0.9_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)
  test_hyper <- nrow(localFileRes)

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates" =c("Covariates1+Covariates2"),
                                   "family_test"=c("gaussian"),
                                   "transformation"="scale",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", markers="DELTAS", areas="GENE")

  fileToRead <- file_path_build(inferenceFolder, "1_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("gaussian"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", markers="DELTAS", areas="GENE")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)



  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("gaussian"),
                                   "transformation"="log10",
                                   "depth_analysis"=2,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential")

  fileToRead <- file_path_build(inferenceFolder, "2_Phenotest_log10_gaussian_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="log",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential")

  fileToRead <- file_path_build(inferenceFolder, "1_Group_log_wilcoxon_", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="quantile_5",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential")

  fileToRead <- file_path_build(inferenceFolder, "1_Group_quantile_5_wilcoxon_", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="quantile_5",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", markers="DELTAS", areas="GENE")
  fileToRead <- file_path_build(inferenceFolder, "3_Group_quantile_5_wilcoxon_", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  })
