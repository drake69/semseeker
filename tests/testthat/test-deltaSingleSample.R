# testthat::test_that("delta_single_sample",{
#
#   
#   tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
#   init_env(tempFolder)
#
#   Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")
#
#   nitem <- 5e4
#   values <- data.frame(Sample_ID=rnorm(nitem, mean=0.5, sd=0.7))
#
#   probe_features <- PROBES[!is.na(PROBES$CHR),]
#   probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
#
#   high_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
#   low_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))
#
#   row.names(values) <- probe_features$PROBE
#   row.names(high_thresholds) <- probe_features$PROBE
#   row.names(low_thresholds) <- probe_features$PROBE
#
#   beta_medians <- high_thresholds - low_thresholds
#
#   delta_single_sample(
#     values = values,
#     high_thresholds = high_thresholds,
#     low_thresholds = low_thresholds,
#     sample_detail = data.frame("Sample_ID"= Sample_ID, "Sample_Group"="Control"),
#     beta_medians = beta_medians,
#     probe_features = probe_features
#   )
#
#   resultFolderData  <-  dir_check_and_create(tempFolder, "Data")
#   outputFolder <- dir_check_and_create(resultFolderData,c("Control","DELTAS_METHYLATION"))
#   fileName <- file_path_build(outputFolder,c(Sample_ID,"DELTAS","METHYLATION"), "bedgraph")
#   expect_true(file.exists(fileName))
#
#
# })
