test_that("analizePopulation", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z]"),sep="")
  init_env(tempFolder)

  nitem <- 5e4
  nsamples <- 20
  methylationData <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylationData <- as.data.frame(matrix(methylationData,nitem,nsamples))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]

  betaSuperiorThresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  betaInferiorThresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(betaSuperiorThresholds) <- probeFeatures$PROBE
  row.names(betaInferiorThresholds) <- probeFeatures$PROBE
  row.names(methylationData) <- probeFeatures$PROBE

  Sample_ID <- stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analizePopulation(methylationData=methylationData,
                    slidingWindowSize = 11,
                    betaSuperiorThresholds = betaSuperiorThresholds,
                    betaInferiorThresholds = betaInferiorThresholds,
                    sampleSheet = mySampleSheet,
                    betaMedians = betaSuperiorThresholds - betaInferiorThresholds,
                    bonferroniThreshold = 0.01,
                    probeFeatures = probeFeatures
                    )

  expect_true(nrow(sp)==nrow(mySampleSheet))
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  # unlink(tempFolder, recursive = TRUE)
  # outputFolder <- dir_check_and_create(resultFolderData,c("Control","MUTATIONS_HYPO"))
  # fileName <- file_path_build(outputFolder,c("MULTIPLE","MUTATIONS","HYPO"), "bed")
  # expect_equal(2 * 2, 4)
})

