test_that("deltaX_get", {

  tempFolder <- tempFolders[1]
  unlink(tempFolder, recursive = TRUE)
  # message(tempFolder)
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
    bonferroni_threshold = bonferroni_threshold,
    inpute="median", start_fresh=TRUE)

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.05

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  keys <- ssEnv$keys_markers_figures_default
  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )
  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Case",],keys)

  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )

  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Control",],keys)

  ss <- semseeker:::deltaX_get(mySampleSheet[mySampleSheet$Sample_Group != "Reference",])

  # verify all DELTAX (except DELTAS and DELTAR ) are coherent with MUTATIONS
  keys <- subset(ssEnv$keys_areas_subareas_markers_figures, AREA == "POSITION")
  keys <- subset(keys, MARKER != "DELTAS")
  keys <- subset(keys, MARKER != "DELTAR")
  keys <- subset(keys, MARKER != "LESIONS")
  keys <- subset(keys, MARKER != "SIGNAL")


  for (k in 1:nrow(keys))
  {
    # k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character(key$AREA)
    subarea <- as.character(key$SUBAREA)

    mutations_pivot_file_name <- semseeker:::pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name))
      mutations_pivot <- polars::pl$read_parquet(mutations_pivot_file_name)$to_data_frame()
    else
      next

    pivot_file_name <- semseeker:::pivot_file_name_parquet(marker,figure,area,subarea)
    testthat::expect_true(file.exists(pivot_file_name))
    if(file.exists(pivot_file_name))
      pivot <- polars::pl$read_parquet(pivot_file_name)$to_data_frame()
    else
      pivot <- NULL


    pivot <- pivot[,-c(1:3)]
    mutations_pivot <- mutations_pivot[,-c(1:3)]
    testthat::expect_true(nrow(pivot)<nprobes)
    testthat::expect_true(nrow(pivot)>0)

    # sort samplesheet by Sample_ID
    mySampleSheet <- mySampleSheet[order(mySampleSheet$Sample_ID),]
    # sort colnames by Sample_ID
    pivot <- pivot[,order(colnames(pivot))]
    mutations_pivot <- mutations_pivot[,order(colnames(mutations_pivot))]

    # all rows have at least one zeo

    testthat::expect_true(all(colnames(pivot) %in% (mySampleSheet$Sample_ID)))
    # wit few positions some sample mutation could be missing
    # testthat::expect_true(all(colnames(pivot) == unique(mySampleSheet$Sample_ID)))
    testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

    testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
    testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

    pivot[is.na(pivot)] <- 0
    mutations_pivot[is.na(mutations_pivot)] <- 0

    pivot[pivot > 0] <- 1
    mutations_pivot[mutations_pivot > 0] <- 1

    if (!all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
    {
      print(marker)
      print(figure)
    }
    testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))

    if(marker!="MUTATIONS")
      for(c in 1:(nrow(mySampleSheet)))
      {
        sample_id <- mySampleSheet[c,"Sample_ID"]
        sample_group <- mySampleSheet[c,"Sample_Group"]
        mutation_bed_file_name <- semseeker:::bed_file_name(sample_id,sample_group,"MUTATIONS",figure)
        if(file.exists(mutation_bed_file_name))
        {
          bed_file_name <- semseeker:::bed_file_name(sample_id,sample_group,marker,figure)
          testthat::expect_true(file.exists(bed_file_name))
        }
      }

  }

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

