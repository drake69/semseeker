test_that("readMultipleBed", {


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

  Sample_ID <- stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
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

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probesPrefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"
  populationName <- unique(Sample_Group)

  probeFeatures <- get(paste0(probesPrefix,"Whole",sep=""))

  res <-readMultipleBed ("MUTATIONS", "BOTH", probeFeatures, columnLabel, populationName, groupingColumnLabel)

  expect_true(nrow(res)>0)

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)

}
)
