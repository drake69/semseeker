test_that("association_analysis", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  nitem <- 6e5
  nsamples <- 30


  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  mySampleSheet$Phenotest <- rnorm(nsamples, mean= 1000, sd= 567)
  mySampleSheet$Group <- c(rep(TRUE,nsamples/2), rep(FALSE,nsamples/2))
  mySampleSheet$Covariates1 <- rnorm(nsamples, mean= 567, sd= 1000)
  mySampleSheet$Covariates2 <- rnorm(nsamples, mean= 67, sd= 100)

  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,parallel_strategy="multisession")
  inferenceFolder <- file.path(tempFolder,"Inference")


  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"="",
                                   "family_test"=c("gaussian"),
                                   "transformation"="scale",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "1_Phenotest_scale_gaussian_test_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)


  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_1000_15000"),
                                   "transformation"="scale",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "1_Phenotest_scale_quantreg_0.5_1000_15000_test_corrected_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("quantreg_0.5_1000_15000"),
                                   "transformation"="scale",
                                   "depth_analysis"=3,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_1000_15000_test_corrected_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Phenotest",
                                   "covariates"=c("Covariates1+Covariates2"),
                                   "family_test"=c("gaussian"),
                                   "transformation"="log10",
                                   "depth_analysis"=2,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "2_Phenotest_log10_gaussian_test_corrected_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="log",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "1_Group_log_wilcoxon_test_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)

  inference_details <- expand.grid("independent_variable"= "Group",
                                   "covariates"="",
                                   "family_test"=c("wilcoxon"),
                                   "transformation"="quantile_5",
                                   "depth_analysis"=1,
                                   "filter_p_value" = FALSE)


  # inference_details,result_folder, maxResources, parallel_strategy
  association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy="multisession")

  fileToRead <- file_path_build(inferenceFolder, "1_Group_quantile_5_wilcoxon_test_result", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  expect_true(nrow(localFileRes)>0)
})

