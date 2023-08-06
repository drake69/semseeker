test_that("create_heatmap", {


  # figures <- c( "BOTH")
  # markers <- c("DELTAS","DELTAR","DELTARQ")
  # areas <- c("GENE")

  # ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures, markers, areas)
  ssEnv <- semseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  # sp <- semseeker:::analyze_population(methylation_data=methylation_data,
  #
  #   sliding_window_size = sliding_window_size,
  #   beta_thresholds = beta_thresholds,
  #   sample_sheet = mySampleSheet,
  #   bonferroni_threshold = bonferroni_threshold,
  #   probe_features = probe_features
  # )

  batch_id <- 1
  sp <- semseeker:::analyze_batch(methylation_data =  methylation_data,
    sample_sheet =  mySampleSheet,
    sliding_window_size = sliding_window_size,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    batch_id = batch_id
  )

  semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)
  sp$Sample_Group <- mySampleSheet$Sample_Group


  sample_groups <- c("Control","Case")
  figures <- c("BOTH")
  markers <- c("DELTAS","DELTAR","DELTARQ","MUTATIONS")

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "PROBE"
  subGroupLabel= "AREA"

  semseeker:::create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  subGroups <- c("CHR")
  probes_prefix = "PROBES_CHR_"
  mainGroupLabel =  "CHR"
  subGroupLabel= "AREA"

  semseeker:::create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  chrBed <- semseeker:::annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
  semseeker:::create_heatmap( inputBedDataFrame =  chrBed,markers = markers, file_prefix = "CHR", groupColumnLabels = c("CHR"))

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/CHR/Control_Vs_Case_CHR_MUTATIONS_BOTH.png")))

  ####################################################################################

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  semseeker:::create_heatmap( inputBedDataFrame = final_bed,markers = markers, file_prefix = "GENE_AREA", groupColumnLabels = c("GENE"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_MUTATIONS_BOTH.png")))

  ####################################################################################

  figures <- c("BOTH")
  markers <- c("DELTAS","DELTAR","DELTARQ","MUTATIONS")

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- semseeker:::annotate_bed (
    sample_groups,
    figures,
    markers,
    subareas,
    probes_prefix,
    area,
    groupingColumnLabel
    )

  semseeker:::create_heatmap( inputBedDataFrame = final_bed,markers = markers, file_prefix = "GENE_AREA", groupColumnLabels = c("AREA"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_DELTAS_BOTH.png")))

  ####################################################################################

  semseeker:::create_heatmap( inputBedDataFrame =  final_bed,markers = markers, file_prefix = "GENE", groupColumnLabels = c("GENE"))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE/Control_Vs_Case_GENE_DELTAS_BOTH.png")))

  ####################################################################################

  # final_bed <- final_bed [1:2,]
  # semseeker:::create_heatmap(inputBedDataFrame = final_bed,markers = markers, file_prefix = "GENE_AREA", groupColumnIDs = c(3))
  #
  # final_bed <- NULL
  # semseeker:::create_heatmap(inputBedDataFrame = final_bed,markers = markers, file_prefix = "GENE_AREA", groupColumnIDs = c(3))

  # unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()
})
