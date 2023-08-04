# test_that("analyze_single_sample_both", {
#
#   library(stringi)
# tmp <- tempdir()
# tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
#   ssEnv <- init_env(tempFolder)
#
#   Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")
#   Sample_Group <- "Control"
#
#   sample_detail <- data.frame("Sample_ID"= Sample_ID,"Sample_Group"= Sample_Group)
#
#   nitem <- 1e4
#   thresholds <- data.frame("thresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
#   values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
# probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
# probe_features <- unique(probe_features)
# probe_features$END <- probe_features$START
# probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
#   row.names(thresholds) <- probe_features$PROBE
#   row.names(values) <- row.names(thresholds)
#
#   mutations <- mutations_get(
#     values = values,
#     figure = "HYPER",
#     thresholds = thresholds,
#     probe_features = probe_features,
#     sampleName = Sample_ID
#   )
#
#   mutations <- subset(mutations, MUTATIONS == 1)[, c("CHR", "START", "END")]
#   folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", "HYPER", sep = "")))
#   dump_sample_as_bed_file(
#     data_to_dump = mutations,
#     fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS","HYPER"),"bed")
#   )
#
#
#   res <- analyze_single_sample(ssEnv =   sample_detail = sample_detail, values = values, thresholds = thresholds, figure = )
#
#   testthat::expect_true(as.numeric(res["MUTATIONS_BOTH"])!=0)
#   testthat::expect_true(ncol(mutations)==3 )
#
#   # da spiegare che test case è quello sotto
#   # testthat::expect_true(mutations[1,2] == mutations[1,3]  )
#
#   fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS","HYPER"),"bed")
#
#   fileData <- read.table(fileName, sep="\t", header = FALSE)
#   testthat::expect_true(nrow(mutations)==nrow(fileData))
#   testthat::expect_true(ncol(fileData)==3 )
#
#   # doParallel::stopImplicitCluster()
#   # parallel::stopCluster(computationCluster)
#
# })
