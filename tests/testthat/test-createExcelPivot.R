test_that("create_excel_pivot", {

  library(stringi)
  tempFolder <- paste("/tmp/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  envir <- init_env(tempFolder, parallel_strategy = "multisession")

  nitem <- 5e3
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

  sp <- analize_population(envir = envir,
                          methylation_data=methylation_data,
                          sliding_window_size = 11,
                          beta_superior_thresholds = beta_superior_thresholds,
                          beta_inferior_thresholds = beta_inferior_thresholds,
                          sample_sheet = sample_sheet,
                          beta_medians = beta_superior_thresholds - beta_inferior_thresholds,
                          bonferroni_threshold = 0.01,
                          probe_features = probe_features
  )

  create_multiple_bed(envir, sample_sheet = sample_sheet, sp)

  populations <- c("Control","Case")

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS","DELTAS")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"

  # create and read
  final_bed <- annotate_bed (
    envir,
    populations ,
    figures ,
    anomalies ,
    groups ,
    probes_prefix ,
    columnLabel ,
    groupingColumnLabel)

  populations <- c("Control","Case")
  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS","DELTAS")

  subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"
  create_excel_pivot (envir=envir, populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  expect_true(file.exists(file.path(envir$result_folderData,"Pivots/GENE.xlsx")))

  subGroups <- c("")
  probes_prefix = "PROBES"
  mainGroupLabel =  "CHR"
  subGroupLabel="GROUP"
  create_excel_pivot (envir=envir, populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

  #TODO: test incremental pivot

#
#   probes_prefix <- "PROBES_Island_"
#   subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
#   mainGroupLabel <- "ISLAND"
#   subGroupLabel <- "RELATION_TO_CPGISLAND"
#   create_excel_pivot ( populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
#
#   subGroups <- c("DMR")
#   probes_prefix = "PROBES_DMR_"
#   mainGroupLabel =  "DMR"
#   subGroupLabel="GROUP"
#   create_excel_pivot (populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)


})
