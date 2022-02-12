test_that("analyseSingleSampleBoth", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  init_env(tempFolder)
  Sample_ID <- stri_rand_strings(1, 7, pattern = "[A-Za-z]")
  Sample_Group <- "Control"

  sampleDetail <- data.frame("Sample_ID"= Sample_ID,"Sample_Group"= Sample_Group)

  nitem <- 8e3
  tresholds <- data.frame("tresholds"= rnorm(nitem, mean=0.5, sd= 0.5))
  values <- data.frame(Sample_ID=rnorm(nitem, mean=0.2, sd=0.5))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  row.names(tresholds) <- probeFeatures$PROBE
  row.names(values) <- row.names(tresholds)

  mutations <- getMutations(
    values = values,
    figure = "HYPER",
    thresholds = tresholds,
    probeFeatures = probeFeatures,
    sampleName = Sample_ID
  )

  mutations <- subset(mutations, MUTATIONS == 1)[, c("CHR", "START", "END")]
  folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("MUTATIONS","_", "HYPER", sep = "")))
  dumpSampleAsBedFile(
    dataToDump = mutations,
    fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"MUTATIONS","HYPER"),"bed")
  )


  res <- analyzeSingleSampleBoth(sampleDetail = sampleDetail)

  expect_true(as.numeric(res["MUTATIONS_BOTH"])!=0)
  expect_true(ncol(mutations)==3 )
  expect_true(mutations[1,2] == mutations[1,3]  )

  fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"MUTATIONS","HYPER"),"bed")

  fileData <- read.table(fileName, sep="\t", header = FALSE)
  expect_true(nrow(mutations)==nrow(fileData))
  expect_true(ncol(fileData)==3 )

})
