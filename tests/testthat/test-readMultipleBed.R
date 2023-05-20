test_that("semseeker:::read_multiple_bed", {


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

  sp$Sample_Group <- mySampleSheet$Sample_Group
  semseeker:::create_multiple_bed(mySampleSheet)

  figures <- c("HYPO", "HYPER", "BOTH")
  anomalies <- c("MUTATIONS","LESIONS","DELTAS","DELTAR","DELTAQ","DELTARQ")

  groups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
  probes_prefix = "PROBES_Gene_"
  columnLabel =  "GENE"
  groupingColumnLabel="GROUP"
  populationName <- "Case"

  probe_features <- probes_get(probes_prefix,"Whole")


  res <-semseeker:::read_multiple_bed ( "MUTATIONS", "BOTH", probe_features, columnLabel, populationName, groupingColumnLabel)
  testthat::expect_true(nrow(res)>0)

  res <-semseeker:::read_multiple_bed ( "DELTAS", "BOTH", probe_features, columnLabel, populationName, groupingColumnLabel)
  testthat::expect_true(nrow(res)>0)
  ####################################################################################

  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPO", probe_features, columnLabel, populationName, groupingColumnLabel)
  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPER", probe_features, columnLabel, populationName, groupingColumnLabel)
  # testthat::expect_true(nrow(res)>0)
  semseeker:::close_env()
}
)
