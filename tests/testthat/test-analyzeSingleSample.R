testthat::test_that("analyzeSingleSample",{

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <-  init_env(tempFolder)
  Sample_ID <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z]")
  Sample_Group <- "Control"

  nitem <- 8e5

  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  row.names(tresholds) <- probeFeatures$PROBE
  row.names(values) <- row.names(tresholds)

  sp <- analyzeSingleSample(envir = envir, values = values,
                      slidingWindowSize = 11,
                      thresholds = tresholds,
                      figure = "HYPO",
                      sampleDetail = data.frame("Sample_ID"=Sample_ID, "Sample_Group"=Sample_Group) ,
                      bonferroniThreshold = 0.05,
                      probeFeatures = probeFeatures)


  outputFolder <- dir_check_and_create(envir$resultFolderData,c("Control","MUTATIONS_HYPO"))
  fileName <- file_path_build(outputFolder,c(Sample_ID,"MUTATIONS","HYPO"), "bed")
  expect_true(file.exists(fileName))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)

})
