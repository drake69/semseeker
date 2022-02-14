test_that("analizePopulation", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
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

  Sample_ID <- stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sampleSheet <- data.frame(Sample_Group, Sample_ID)
  betaMedians <- betaSuperiorThresholds + betaInferiorThresholds / 2
  slidingWindowSize <- 11
  bonferroniThreshold <- 0.01

  sp <- analizePopulation(methylationData=methylationData,
                    slidingWindowSize = slidingWindowSize,
                    betaSuperiorThresholds = betaSuperiorThresholds,
                    betaInferiorThresholds = betaInferiorThresholds,
                    sampleSheet = sampleSheet,
                    betaMedians = betaMedians,
                    bonferroniThreshold = bonferroniThreshold,
                    probeFeatures = probeFeatures
                    )

  message(nrow(sp))
  message(nrow(sampleSheet))
  expect_true(nrow(sp)==nrow(sampleSheet))
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)
  # unlink(tempFolder, recursive = TRUE)
  # outputFolder <- dir_check_and_create(resultFolderData,c("Control","MUTATIONS_HYPO"))
  # fileName <- file_path_build(outputFolder,c("MULTIPLE","MUTATIONS","HYPO"), "bed")
  # expect_equal(2 * 2, 4)
})

