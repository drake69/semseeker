test_that("association_analysis", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  nitem <- 1e3
  nsamples <- 30

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))

  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))
  row.names(methylation_data) <- probe_features$PROBE[1:nitem]

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  mySampleSheet$Phenotest <- rnorm(nsamples, mean= 1000, sd= 567)
  mySampleSheet$Group <- c(rep(TRUE,nsamples/2), rep(FALSE,nsamples/2))
  mySampleSheet$Covariates1 <- rnorm(nsamples, mean= 567, sd= 1000)
  mySampleSheet$Covariates2 <- rnorm(nsamples, mean= 67, sd= 100)

  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,
    parallel_strategy="sequential", figures="BOTH", anomalies="DELTAS", metaareas="PROBE")
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
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="PROBE")

  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("quantreg_0.5"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential",
    figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_100_1000_0.99"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="PROBE")

  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("quantreg_0.5_100_1000_0.99_np"),
    "transformation"="scale",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="PROBE")


  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,parallel_strategy="sequential", figures="BOTH", anomalies="DELTAS", metaareas="CHR")
  inferenceFolder <- file.path(tempFolder,"Inference")


  #todo: test incremental association analysis

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_1000_15000_0.9"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="BOTH", anomalies=c("DELTAS","DELTAQ"), metaareas="CHR")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_1000_15000_0.9_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  test_both <- nrow(localFileRes)
  testthat::expect_true(nrow(localFileRes)>0)

  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", anomalies="DELTAS", metaareas="CHR")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_1000_15000_0.9_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)
  test_hyper <- nrow(localFileRes)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="log",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", metaareas="GENE")

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
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", metaareas="GENE")

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
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", anomalies="DELTAS", metaareas="GENE")
  fileToRead <- file_path_build(inferenceFolder, "3_Group_quantile_5_wilcoxon_", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates" =c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="scale",
  #   "depth_analysis"=1,
  #   "filter_p_value" = FALSE)
  #
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", anomalies="DELTAS", metaareas="GENE")
  #
  # fileToRead <- file_path_build(inferenceFolder, "1_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="scale",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  #
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", figures="HYPER", anomalies="DELTAS", metaareas="GENE")
  #
  # fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  #
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="log10",
  #   "depth_analysis"=2,
  #   "filter_p_value" = FALSE)
  #
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="sequential", metaareas="GENE")
  #
  # fileToRead <- file_path_build(inferenceFolder, "2_Phenotest_log10_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  })
