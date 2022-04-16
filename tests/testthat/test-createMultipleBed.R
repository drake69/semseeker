test_that("createMultipleBed", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

  nitem <- 5e4
  nsamples <- 5
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))
  probe_features <- PROBES[!is.na(PROBES$CHR),]
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

  createMultipleBed(envir, sample_sheet)

  tempresultFolder <-dir_check_and_create(envir$resultFolderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresultFolder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "bed")
  localFileRes <- read.table(fileToRead, sep="\t")

  expect_true(nrow(localFileRes)>0)

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)


})
