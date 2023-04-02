test_that("population_check", {


  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

  nitem <- 1e3
  nsamples <- 30
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)


  # all fine
  expect_true(is.null(population_check(mySampleSheet, methylation_data, envir)))

  #Sample_ID has NA
  mySampleSheet$Sample_ID[1] <- NA
  expect_true(!is.null(population_check(mySampleSheet, methylation_data, envir)))

  #Sample_Group has NA
  mySampleSheet$Sample_Group[1] <- NA
  expect_true(!is.null(population_check(mySampleSheet, methylation_data, envir)))

  #Lost Sample_Group Values
  mySampleSheet$Sample_Group <- NA
  expect_true(!is.null(population_check(mySampleSheet, methylation_data, envir)))

  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(computationCluster)

})

