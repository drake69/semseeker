test_that("semeeker", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "multisession")


  nitem <- 5e3
  nsamples <- 21

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data_1 <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data_1 <- as.data.frame(matrix(methylation_data_1,nitem,nsamples))
  row.names(methylation_data_1) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data_1) <- Sample_ID

  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet_1 <- data.frame(Sample_Group, Sample_ID)

  #########

  nitem <- 5e3
  nsamples <- 21

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data_2 <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data_2 <- as.data.frame(matrix(methylation_data_2,nitem,nsamples))
  row.names(methylation_data_2) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data_2) <- Sample_ID

  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet_2 <- data.frame(Sample_Group, Sample_ID)

  ##########


  nitem <- 5e3
  nsamples <- 21

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START

  nitem <- min(nitem, nrow(probe_features))
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  methylation_data_3 <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data_3 <- as.data.frame(matrix(methylation_data_3,nitem,nsamples))
  row.names(methylation_data_3) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data_3) <- Sample_ID

  Sample_Group <- c(rep("Control",nsamples/3),rep("Reference",nsamples/3),rep("Case",nsamples/3))
  mySampleSheet_3 <- data.frame(Sample_Group, Sample_ID)


  row.names(methylation_data_2) <- row.names(methylation_data_1)
  row.names(methylation_data_3) <- row.names(methylation_data_1)

  mySampleSheet <- list(mySampleSheet_1, mySampleSheet_2, mySampleSheet_3)
  methylation_data <- list(methylation_data_1, methylation_data_2, methylation_data_3)

  semseeker( sample_sheet =  mySampleSheet,methylation_data =  methylation_data, result_folder = tempFolder,
             parallel_strategy = "multisession", anomalies="DELTAQ", metaareas="GENE", figures="BOTH")

  # batch_correlation_check(envir)
  tempresult_folder <- file.path(tempFolder,"Data","Control","MUTATIONS_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "MUTATIONS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  # localFileRes <- read.table(fileToRead, sep="\t")

  testthat::expect_true(nrow(localFileRes_both)>0)

  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAS_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "DELTAS" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  # localFileRes <- read.table(fileToRead, sep="\t")

  testthat::expect_true(nrow(localFileRes_both)>0)


  # test deltaq creation
  tempresult_folder <- file.path(tempFolder,"Data","Control","DELTAQ_BOTH")
  fileToRead <- file_path_build(tempresult_folder, c("MULTIPLE", "DELTAQ" ,"BOTH" ), "fst")
  localFileRes_both <- fst::read_fst(fileToRead)
  testthat::expect_true(sum(is.na(localFileRes_both$VALUE))==0)

})

