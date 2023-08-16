test_that("beta",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  sample_detail <- mySampleSheet[1,]

  res <- semseeker:::beta_single_sample(values = methylation_data[,sample_detail$Sample_ID], sample_detail = sample_detail, probe_features = probe_features)

  folder_to_save <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c(sample_detail$Sample_Group,paste0("BETA","_", "MEAN", sep = "")))
  beta_file <- file_path_build(folder_to_save,c(sample_detail$Sample_ID,"BETA","MEAN"),"bedgraph")
  testthat::expect_true(file.exists(beta_file))

  beta_file <- read.table(beta_file, sep="\t")
  testthat::expect_true(nrow(beta_file)==nrow(methylation_data))

})
