test_that("populationCheck", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stri_rand_strings(1, 7, pattern = "[A-Za-z]"),sep="")
  init_env(tempFolder)

  nitem <- 5
  nsamples <- 30
  methylationData <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylationData <- as.data.frame(matrix(methylationData,nitem,nsamples))

  Sample_ID <- stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylationData) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)


  # all fine
  expect_true(is.null(populationCheck(mySampleSheet, methylationData)))

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  expect_true(!is.null(populationCheck(mySampleSheet, methylationData)))

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  expect_true(!is.null(populationCheck(mySampleSheet, methylationData)))

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  expect_true(!is.null(populationCheck(mySampleSheet, methylationData)))
})

