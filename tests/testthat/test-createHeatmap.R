test_that("createHeatmap", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

  nitem <- 5e5
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

  sp <- analizePopulation(methylationData=methylationData,
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
    populations ,
    figures ,
    anomalies ,
    groups ,
    probesPrefix ,
    columnLabel ,
    groupingColumnLabel)

  createHeatmap(inputBedDataFrame = finalBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnIDs = c(3))
  expect_true(file.exists(file.path(envir$resultFolderChart,"/GENE_AREA/Control_GENE_AREA_MUTATIONS.png")))

  # finalBed <- finalBed [1:2,]
  # createHeatmap(inputBedDataFrame = finalBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnIDs = c(3))
  #
  # finalBed <- NULL
  # createHeatmap(inputBedDataFrame = finalBed,anomalies = anomalies, groupLabel = "GENE_AREA", groupColumnIDs = c(3))

})
