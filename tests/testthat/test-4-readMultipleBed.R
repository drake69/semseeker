test_that("semseeker:::read_multiple_bed", {


  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(methylation_data)

  ####################################################################################
  sp <- semseeker:::analyze_population(methylation_data=methylation_data,

    sliding_window_size = sliding_window_size,
    beta_thresholds = beta_thresholds,
    sample_sheet = mySampleSheet,
    bonferroni_threshold = bonferroni_threshold,
    probe_features = probe_features

  )

  sample_group <- unique(mySampleSheet$Sample_Group)[1]
  semseeker:::create_multiple_bed(mySampleSheet)



  res <-semseeker:::read_multiple_bed (marker =  "MUTATIONS",figure =  "BOTH", sample_group = sample_group)
  testthat::expect_true(nrow(res)>0)

  res <-semseeker:::read_multiple_bed (marker =  "DELTAS",figure =  "BOTH", sample_group = sample_group)
  testthat::expect_true(nrow(res)>0)
  ####################################################################################

  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPO", probe_features, area, sample_group, groupingColumnLabel)
  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPER", probe_features, area, sample_group, groupingColumnLabel)
  # testthat::expect_true(nrow(res)>0)
  semseeker:::close_env()
}
)
