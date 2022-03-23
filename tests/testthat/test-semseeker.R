test_that("semeeker", {

  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

  nitem <- 5^4
  nsamples <- 30


  methylationData <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylationData <- as.data.frame(matrix(methylationData,nitem,nsamples))
  probeFeatures <- PROBES[!is.na(PROBES$CHR),]
  probeFeatures <- probeFeatures[probeFeatures$PROBE %in% sample(x=probeFeatures[,"PROBE"] , size=nitem),]
  row.names(methylationData) <- probeFeatures$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  semseeker( sampleSheet =  mySampleSheet,methylationData =  methylationData, resultFolder = tempFolder)

  tempresultFolder <-dir_check_and_create(envir$resultFolderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresultFolder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "bed")
  localFileRes <- read.table(fileToRead, sep="\t")

  expect_true(nrow(localFileRes)>0)

})

