test_that("semeeker", {

  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder)

  nitem <- 5^4
  nsamples <- 30


  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))
  probe_features <- PROBES[!is.na(PROBES$CHR),]
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder)

  tempresult_folder <-dir_check_and_create(envir$result_folderData,c("Control","MUTATIONS_BOTH"))
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "bed")
  localFileRes <- read.table(fileToRead, sep="\t")

  expect_true(nrow(localFileRes)>0)

})

