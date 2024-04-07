test_that("semseeker:::read_multiple_bed", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  semseeker:::create_multiple_bed(mySampleSheet)
  dq <- deltaq_get(mySampleSheet)
  drq <- deltarq_get(mySampleSheet)

  markers <-c("MUTATIONS","DELTAS","DELTAQ","DELTARQ","DELTAR")
  figures <- c("HYPO","HYPER","BOTH")
  results <- data.frame("MARKER","FIGURE","SAMPLE_GROUP","NROWS")
  results <- results[0,]

  for (sg in length(unique(mySampleSheet$Sample_Group)))
  {
    # sg <-1
    sample_group <- unique(mySampleSheet$Sample_Group)[sg]
    for (marker in markers)
    {
      # marker <- "DELTARQ"
      for (figure in figures)
      {
        # figure <- "HYPO"
        # message(sample_group,figure,marker)
        res <-semseeker:::read_multiple_bed (marker =  marker,figure =  figure, sample_group = sample_group)
        testthat::expect_true(nrow(res)>0)
        results <- rbind(results, data.frame("MARKER"=marker,"FIGURE"=figure,"SAMPLE_GROUP"=sample_group,"NROWS"=nrow(res)))
      }
    }
  }


  # create combination of markers
  comb <- combn(markers,2)
  for (i in 1:ncol(comb)) {
    # i <- 1
    for (figure in figures)
    {
      # figure <- "HYPO"
      for (sample_group in unique(mySampleSheet$Sample_Group))
      {
        res2 <- semseeker:::read_multiple_bed (comb[2,i], figure, sample_group)
        res1 <- semseeker:::read_multiple_bed (marker =  comb[1,i], figure =  figure, sample_group =  sample_group)
        if(comb[2,i]=="MUTATIONS")
        {
          res2 <- nrow(res2)
        } else
        {
          res2 <- sum(res2$VALUE)
        }
        if(comb[1,i]=="MUTATIONS")
        {
          res1 <- nrow(res1)
        } else
        {
          res1 <- sum(res1$VALUE)
        }
        # message(comb[2,i]," ",comb[1,i]," ", sample_group, " ",figure, " ", res1!=res2)
        testthat::expect_true(res1!=res2)
      }
    }
  }

  ####################################################################################

  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPO", probe_features, area, sample_group, groupingColumnLabel)
  # res <-semseeker:::read_multiple_bed ( "DELTAS", "HYPER", probe_features, area, sample_group, groupingColumnLabel)
  # testthat::expect_true(nrow(res)>0)
  semseeker:::close_env()
}
  )
