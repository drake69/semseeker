test_that("create_excel_pivot", {

  ssEnv <- init_env(tempFolder, parallel_strategy = "sequential")

  nitem <- 1e4
  nsamples <- 10
  methylation_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  methylation_data <- as.data.frame(matrix(methylation_data,nitem,nsamples))

  probe_features <- PROBES_Gene_Whole[!is.na(PROBES_Gene_Whole$START),c("CHR","START","PROBE")]
  probe_features <- unique(probe_features)
  probe_features$END <- probe_features$START
  probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]

  beta_superior_thresholds <- data.frame(rnorm(nitem, mean = 1, sd=0.2))
  beta_inferior_thresholds <- data.frame(rnorm(nitem, mean=0.2, sd=0.2))

  row.names(beta_superior_thresholds) <- probe_features$PROBE
  row.names(beta_inferior_thresholds) <- probe_features$PROBE
  row.names(methylation_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 7, pattern = "[A-Za-z]")
  colnames(methylation_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/2), rep("Case",nsamples/2))
  sample_sheet <- data.frame(Sample_Group, Sample_ID)

  sp <- analyze_population(methylation_data=methylation_data,
    sliding_window_size = 11,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features,
    bonferroni_threshold = 0.01,
  )
  sp$Sample_Group <- sample_sheet$Sample_Group

  create_multiple_bed( sample_sheet = sample_sheet)

  sample_groups <- c("Control","Case")

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("MUTATIONS","LESIONS","DELTAS")

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  area =  "GENE"
  groupingColumnLabel="AREA"

  # create and read
  final_bed <- annotate_bed (

    sample_groups ,
    figures ,
    markers ,
    subareas ,
    probes_prefix ,
    area ,
    groupingColumnLabel)

  sample_groups <- c("Control","Case")
  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("MUTATIONS","LESIONS","DELTAS")

  subGroups <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  probes_prefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="AREA"
  create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/GENE.xlsx")))

  subGroups <- c("")
  probes_prefix = "PROBES_CHR_"
  mainGroupLabel =  "CHR"
  subGroupLabel="AREA"
  create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "PROBE"
  subGroupLabel="AREA"
  create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  #TODO: test incremental pivot

#
#   probes_prefix <- "PROBES_Island_"
#   subGroups <- c("N_SHORE","S_SHORE","N_Shelf","S_SHELF", "WHOLE")
#   mainGroupLabel <- "ISLAND"
#   subGroupLabel <- "RELATION_TO_CPGISLAND"
#   create_excel_pivot ( sample_groups, figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
#
#   subGroups <- c("DMR")
#   probes_prefix = "PROBES_DMR_"
#   mainGroupLabel =  "DMR"
#   subGroupLabel="AREA"
#   create_excel_pivot (sample_groups, figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)


})
