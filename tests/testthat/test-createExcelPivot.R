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

  populations <- c("Control","Case")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("DELTAR","DELTARQ","MUTATIONS","LESIONS","DELTAS","DELTAQ","BETA")
  # figures <- c("HYPO")
  # anomalies <- c("DELTAR")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  # groups <- c("Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- semseeker:::annotate_bed (
    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  populations <- c("Control","Case")
  # figures <- c("HYPO", "HYPER", "BOTH")
  # anomalies <- c("DELTAR","DELTARQ","MUTATIONS","LESIONS","DELTAS","DELTAQ")

  # subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  subGroups <- groups
  probes_prefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"
  semseeker:::create_excel_pivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/GENE.xlsx")))

  ####################################################################################

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "PROBE"
  subGroupLabel="GROUP"
  semseeker:::create_excel_pivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/PROBE.xlsx")))

  ####################################################################################

    # subGroups <- c("")
  # probes_prefix = "PROBES_CHR_"
  # mainGroupLabel =  "CHR"
  # subGroupLabel="GROUP"
  # semseeker:::create_excel_pivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
  # testthat::expect_true(file.exists(file.path(ssEnv$result_folderData,"Pivots/CHR.xlsx")))

  #TODO: test incremental pivot

#
#   probes_prefix <- "PROBES_Island_"
#   subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
#   mainGroupLabel <- "ISLAND"
#   subGroupLabel <- "RELATION_TO_CPGISLAND"
#   semseeker:::create_excel_pivot ( populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
#
#   subGroups <- c("DMR")
#   probes_prefix = "PROBES_DMR_"
#   mainGroupLabel =  "DMR"
#   subGroupLabel="GROUP"
#   semseeker:::create_excel_pivot (populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)

  semseeker:::close_env()
})
