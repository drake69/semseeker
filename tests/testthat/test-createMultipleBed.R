test_that("create_multiple_bed", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e3
  nsamples <- 5
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(beta_superior_thresholds) <- probe_features$PROBE
  row.names(beta_inferior_thresholds) <- probe_features$PROBE
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sample_sheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analize_population(envir,
                           methylation_data=methylation_data,
                           sliding_window_size = 11,
                           beta_superior_thresholds = beta_superior_thresholds,
                           beta_inferior_thresholds = beta_inferior_thresholds,
                           sample_sheet = sample_sheet,
                           beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                           bonferroni_threshold = 0.01,
                           probe_features = probe_features
  )

  sp$Sample_Group <- sample_sheet$Sample_Group
  create_multiple_bed(envir, sample_sheet)

  tempresult_folder <-dir_check_and_create(envir$result_folderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  # localFileRes_both <- read.table(fileToRead, sep="\t")
  localFileRes_both <- fst::read_fst(fileToRead)

  expect_true(nrow(localFileRes_both)>0)

  # tempresult_folder <-dir_check_and_create(envir$result_folderData,c("Control","MUTATIONS_HYPO"))
  # fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"HYPO" ), "fst")
  # localFileRes_hypo <- read.table(fileToRead, sep="\t")
  #
  # tempresult_folder <-dir_check_and_create(envir$result_folderData,c("Control","MUTATIONS_HYPER"))
  # fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"HYPER" ), "fst")
  # localFileRes_hyper <- read.table(fileToRead, sep="\t")

  # expect_true(nrow(localFileRes_hyper)>0 | nrow(localFileRes_hypo)>0 | nrow(localFileRes_both)>0)

  close_env(envir)
})
