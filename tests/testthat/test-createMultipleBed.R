test_that("createMultipleBed", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder)

  nitem <- 5e4
  nsamples <- 5
  methylationData <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylationData <- as.data.frame(matrix(methylationData,nitem,nsamples))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]

  betaSuperiorThresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  betaInferiorThresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(betaSuperiorThresholds) <- probeFeatures$PROBE
  row.names(betaInferiorThresholds) <- probeFeatures$PROBE
  row.names(methylationData) <- probeFeatures$PROBE

  Sample_ID <- stri_rand_strings(nsamples, 7, pattern = "[A-Za-z0-9]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sampleSheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analizePopulation(methylationData=methylationData,
                          slidingWindowSize = 11,
                          betaSuperiorThresholds = betaSuperiorThresholds,
                          betaInferiorThresholds = betaInferiorThresholds,
                          sampleSheet = sampleSheet,
                          betaMedians = betaSuperiorThresholds - betaInferiorThresholds,
                          bonferroniThreshold = 0.01,
                          probeFeatures = probeFeatures
  )

  createMultipleBed(sampleSheet)

  tempresultFolder <-dir_check_and_create(resultFolderData,c("Control","MUTATIONS_BOTH"))

  fileToRead <- file_path_build(tempresultFolder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "bed")
  localFileRes <- read.table(fileToRead, sep="\t")

  expect_true(nrow(localFileRes)>0)


})
