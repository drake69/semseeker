test_that("annotateBed", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

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

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sampleSheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analizePopulation(envir,
                          methylationData=methylationData,
                          slidingWindowSize = 11,
                          betaSuperiorThresholds = betaSuperiorThresholds,
                          betaInferiorThresholds = betaInferiorThresholds,
                          sampleSheet = sampleSheet,
                          betaMedians = betaSuperiorThresholds - betaInferiorThresholds,
                          bonferroniThreshold = 0.01,
                          probeFeatures = probeFeatures
  )

  createMultipleBed(envir, sampleSheet = sampleSheet)

  populations <- c("Control")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probesPrefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  finalBed <- annotateBed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probesPrefix ,
    columnLabel ,
    groupingColumnLabel)

  bedFileName <- file_path_build(envir$resultFolderData , c(columnLabel, "ANNOTATED"),"csv")

  # file extsits
  expect_true(file.exists(bedFileName))

  # not empty data set
  expect_true(nrow(finalBed)>0)

  # has the correct header
  expect_true( columnLabel %in% colnames(finalBed))

  #read again  existent
  finalBed <- annotateBed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probesPrefix ,
    columnLabel ,
    groupingColumnLabel)

  expect_true( columnLabel %in% colnames(finalBed))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)


})
