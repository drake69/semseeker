testthat::test_that("analyze_single_sample",{

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <-  init_env(tempFolder, parallel_strategy = "sequential")
  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")
  Sample_Group <- "Control"

  nitem <- 1e3

  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))

  probes <- probes_get("PROBES_Gene_","Whole")
  probe_features <- probes[!is.na(probes$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  row.names(tresholds) <- probe_features$PROBE
  row.names(values) <- row.names(tresholds)

  ####################################################################################

  sp <- analyze_single_sample(envir = envir, values = values,
                      sliding_window_size = 11,
                      thresholds = tresholds,
                      figure = "HYPO",
                      sample_detail = data.frame("Sample_ID"=Sample_ID, "Sample_Group"=Sample_Group) ,
                      bonferroni_threshold = 0.05,
                      probe_features = probe_features)


  outputFolder <- dir_check_and_create(envir$result_folderData,c("Control","MUTATIONS_HYPO"))
  fileName <- file_path_build(outputFolder,c(Sample_ID,"MUTATIONS","HYPO"), "bed")
  expect_true(file.exists(fileName))

  ####################################################################################

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  close_env(envir)
})
