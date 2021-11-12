test_that("readMultipleBed", {


  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder)
  Sample_ID <- stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")
  Sample_Group <- "Control"

  nitem <- 100

  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  row.names(tresholds) <- probeFeatures$PROBE
  row.names(values) <- row.names(tresholds)

  sampleDetail <- data.frame("Sample_ID"=Sample_ID, "Sample_Group"=Sample_Group)

  sp <- analyzeSingleSample(values = values,
                            slidingWindowSize = 11,
                            thresholds = tresholds,
                            figure = "HYPO",
                            sampleDetail = sampleDetail ,
                            bonferroniThreshold = 0.05,
                            probeFeatures = probeFeatures)


  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probesPrefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"
  populationName <- Sample_Group

  probeFeatures <- get(paste0(probesPrefix,"Whole",sep=""))



  res <-readMultipleBed ("MUTATIONS", "HYPO", probeFeatures, columnLabel, populationName, groupingColumnLabel)

  # outputFolder <- dir_check_and_create(resultFolderData,c("Control","MUTATIONS_HYPO"))
  # fileName <- file_path_build(outputFolder,c(Sample_ID,"MUTATIONS","HYPO"), "bed")
  expect_true(nrow(res)>0)
}
)
