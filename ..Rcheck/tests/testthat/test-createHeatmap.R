test_that("create_heatmap", {

  library(stringi)
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  figures <- c( "BOTH")
  anomalies <- c("DELTAS")
  metaareas <- c("GENE")

  ssEnv <- init_env(result_folder =  tempFolder, parallel_strategy = "sequential", maxResources = 90, figures, anomalies, metaareas)

  nitem <- 1e4
  nsamples <- 5

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]


  nitem <- min(nitem, nrow(probe_features))
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(beta_superior_thresholds) <- probe_features$PROBE
  row.names(beta_inferior_thresholds) <- probe_features$PROBE
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- rep("Control",nsamples)
  sample_sheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analize_population(methylation_data=methylation_data,
    sliding_window_size = 11,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features,
    bonferroni_threshold = 0.01,
  )

  create_multiple_bed( sample_sheet = sample_sheet)
  sp$Sample_Group <- sample_sheet$Sample_Group


  populations <- c("Control")
  figures <- c("BOTH")
  anomalies <- c("MUTATIONS")

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "PROBE"
  subGroupLabel= "GROUP"

  create_excel_pivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  subGroups <- c("CHR")
  probes_prefix = "PROBES_CHR_"
  mainGroupLabel =  "CHR"
  subGroupLabel= "GROUP"

  create_excel_pivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  chrBed <- annotate_bed(populations ,figures ,anomalies ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
  create_heatmap( inputBedDataFrame =  chrBed,anomalies = anomalies, file_prefix = "CHR", groupColumnLabels = c("CHR"))

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/CHR/Control_CHR_MUTATIONS_BOTH.png")))

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  create_heatmap( inputBedDataFrame = final_bed,anomalies = anomalies, file_prefix = "GENE_AREA", groupColumnLabels = c("GENE"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_GENE_AREA_MUTATIONS_BOTH.png")))

  figures <- c("BOTH")
  anomalies <- c("DELTAS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (

    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  create_heatmap( inputBedDataFrame = final_bed,anomalies = anomalies, file_prefix = "GENE_AREA", groupColumnLabels = c("GROUP"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_GENE_AREA_DELTAS_BOTH.png")))


  create_heatmap( inputBedDataFrame =  final_bed,anomalies = anomalies, file_prefix = "GENE", groupColumnLabels = c("GENE"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE/Control_GENE_DELTAS_BOTH.png")))


  # final_bed <- final_bed [1:2,]
  # create_heatmap(inputBedDataFrame = final_bed,anomalies = anomalies, file_prefix = "GENE_AREA", groupColumnIDs = c(3))
  #
  # final_bed <- NULL
  # create_heatmap(inputBedDataFrame = final_bed,anomalies = anomalies, file_prefix = "GENE_AREA", groupColumnIDs = c(3))

})
