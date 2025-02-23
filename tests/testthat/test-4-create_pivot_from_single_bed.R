test_that("create_position_pivot_from_single_bed", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  ####################################################################################

  semseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data)
  }
  probe_features <<- semseeker::PROBES[semseeker::PROBES$PROBE %in% rownames(signal_data),]

  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )

  sp <- semseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )

  ssEnv <- semseeker:::get_session_info()

  keys <- expand.grid(
    MARKER = c("DELTAS","DELTAR","MUTATIONS","LESIONS"),
    FIGURE = c("HYPER","HYPO","BOTH")
  )
  area <- "PROBE"
  subarea <- "WHOLE"

  for (k in 1:nrow(keys))
  {

    key <- keys[k,]
    marker <- key$MARKER
    figure <- key$FIGURE

    for ( s in 1:nrow(mySampleSheet))
    {
      sample <- mySampleSheet[s,]
      create_position_pivot_from_single_bed(sample$Sample_Group,sample$Sample_ID, marker,figure,"PROBE","WHOLE")
    }
    reportFolder <- dir_check_and_create(ssEnv$result_folderData,"PivotsNew")
    pivot_subfolder <- dir_check_and_create(reportFolder, marker)
    pivot_file_name <- paste0(marker,"_", figure,"_",area,"_",subarea, sep = "")
    pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_file_name,extension =  ".csv" ,add_gz=TRUE)
    testthat::expect_true(file.exists(pivot_file_name))
    if(file.exists(pivot_file_name))
      pivot <- readr::read_delim(pivot_file_name, show_col_types=FALSE, progress=FALSE)
    else
      pivot <- NULL
    testthat::expect_true(ncol(pivot)==(nsamples+3))
    testthat::expect_true(nrow(pivot)<nprobes)
    testthat::expect_true(nrow(pivot)>0)
  }


  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
