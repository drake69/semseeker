test_that("get-meth_tech", {

  
  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  figures <- c( "BOTH")
  anomalies <- c("DELTAS","DELTAQ")
  metaareas <- c("GENE")

  init_env(result_folder =  tempFolder, parallel_strategy = "sequential", maxResources = 90, figures = "BOTH", anomalies = "DELTAS", metaareas = "GENE")
  # envir <- .pkgglobalenv$ssEnv
  # expect_true(length(envir$result_folderData)>0)

  nitem <- 1e3
  nsamples <- 5

  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probes <- semseeker::PROBES
  probes <- subset(probes, k450 & ! k850)

  probe_features <- probes[!is.na(probes$START),c("CHR","START","PROBE")]
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
  row.names(methylation_data) <- probe_features$PROBE

  ssEnv <- get_meth_tech(methylation_data)
  # expect_true(length(envir$result_folderData)>0)
  expect_true(ssEnv$tech=="k450")
  # expect_true(envir$tech,"k850")

  # expect_equal(.pkgglobalenv$myvar, 42)

})
