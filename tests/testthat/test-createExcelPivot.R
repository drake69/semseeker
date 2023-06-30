test_that("create_excel_pivot", {

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################

  sp <- semseeker:::analize_population(methylation_data=methylation_data,
    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features

  )


  semseeker:::create_multiple_bed( sample_sheet = mySampleSheet)

  ####################################################################################

  sample_groups <- c("Control","Case")

  figures <- c("HYPO", "HYPER", "BOTH")
  markers <- c("DELTAR","DELTARQ","MUTATIONS","LESIONS","DELTAS","DELTAQ","BETA")
  # figures <- c("HYPO")
  # markers <- c("DELTAR")

  subareas <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  # subareas <- c("WHOLE")
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

  sample_groups <- c("Control","Case")
  # figures <- c("HYPO", "HYPER", "BOTH")
  # markers <- c("DELTAR","DELTARQ","MUTATIONS","LESIONS","DELTAS","DELTAQ")

  # subGroups <- c("BODY","TSS1500","5UTR","TSS200","1STEXON","3UTR","EXNBND","WHOLE")
  subGroups <- subareas
  probes_prefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="AREA"
  semseeker:::create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/GENE.xlsx")))

  ####################################################################################

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "PROBE"
  subGroupLabel="AREA"
  semseeker:::create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/PROBE.xlsx")))

  ####################################################################################

    # subGroups <- c("")
  # probes_prefix = "PROBES_CHR_"
  # mainGroupLabel =  "CHR"
  # subGroupLabel="AREA"
  # semseeker:::create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  # testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/CHR.xlsx")))

  #TODO: test incremental pivot

#
#   probes_prefix <- "PROBES_Island_"
#   subGroups <- c("N_SHORE","S_SHORE","N_SHELF","S_SHELF","WHOLE")
#   mainGroupLabel <- "ISLAND"
#   subGroupLabel <- "RELATION_TO_CPGISLAND"
#   semseeker:::create_excel_pivot ( sample_groups, figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
#
#   subGroups <- c("DMR")
#   probes_prefix = "PROBES_DMR_"
#   mainGroupLabel =  "DMR"
#   subGroupLabel="AREA"
#   semseeker:::create_excel_pivot (sample_groups, figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)

  semseeker:::close_env()
})
