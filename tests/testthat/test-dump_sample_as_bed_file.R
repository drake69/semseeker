# test_that("semseeker:::dump_sample_as_bed_file", {
#
#
#   tmp <- tempdir()
#   tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
#   ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)
#
#   sample_detail <- mySampleSheet[1,]
#   beta_values <- as.data.frame(methylation_data[, sample_detail$Sample_ID])
#
#   semseeker:::dump_sample_as_bed_file(
#     data_to_dump = beta_values,
#     fileName = semseeker:::file_path_build(folder_to_save,c(sample_detail$Sample_ID,"BETAS"),"bed")
#   )
#
#   res <-semseeker:::read_multiple_bed ( "MUTATIONS", "BOTH", probe_features, area, sample_group, groupingColumnLabel)
#   testthat::expect_true(nrow(res)>0)
#
#   res <-semseeker:::read_multiple_bed ( "DELTAS", "BOTH", probe_features, area, sample_group, groupingColumnLabel)
#   testthat::expect_true(nrow(res)>0)
#   ####################################################################################
#
#   # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPO", probe_features, area, sample_group, groupingColumnLabel)
#   # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPER", probe_features, area, sample_group, groupingColumnLabel)
#   # testthat::expect_true(nrow(res)>0)
#   semseeker:::close_env()
# }
# )
