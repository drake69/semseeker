#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param sample_sheet dataframe with at least a column Sample_ID to identify samples
#' @param methylation_data matrix of methylation data
#' @param bonferroni_threshold = 0.05 #threshold to define which pValue adjusted to define an epilesion
#' @param maxResources percentage of how many available cores will be used default 90 percent, rounded to the lowest integer
#' @param iqrTimes how many times below the first quartile and over the third quartile the interqauartile is "added" to define the outlier
#' @param parallel_strategy which strategy to use for parallel executio see future vignete: possibile values, none, multisession,sequential, multicore, cluster
#' @param result_folder where the result will be saved
#' @param ... other options to filter elaborations
#'
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#' @importFrom doRNG %dorng%
#'
semseeker <- function(sample_sheet,
                      methylation_data,
                      result_folder,
                      bonferroni_threshold = 0.05,
                      maxResources = 90,
                      iqrTimes = 3 ,
                      parallel_strategy ="multisession",
                      ... ) {


  unlink(result_folder, recursive = TRUE)
  init_env( result_folder= result_folder, maxResources= maxResources, parallel_strategy = parallel_strategy, ...)
  ssEnv <- get_session_info()

  # set digits to 22
  withr::local_options(list(digits = 22))
  sliding_window_size <- 11

  if(is.data.frame(sample_sheet) & is.data.frame(methylation_data))
  {
    sample_sheet <-list(sample_sheet)
    methylation_data <- list(methylation_data)
  } else
  {
    if(is.data.frame(sample_sheet) | is.data.frame(methylation_data))
      stop("both sample_sheet and methylation_data should be data frame!")
    if(length(sample_sheet)!=length(methylation_data))
      stop("both sample_sheet and methylation_data should have been list with the same length!")
  }

  if(length(methylation_data)>1)
  {
    d <- 1
    for(d in 1:length(methylation_data))
    {
      if(exists("probes_to_preserve"))
        probes_to_preserve <- probes_to_preserve[ probes_to_preserve %in% row.names(stats::na.omit(methylation_data[[d]]))]
      else
        probes_to_preserve <- row.names(stats::na.omit(methylation_data[[d]]))
    }
  }
  else
    probes_to_preserve <- row.names(methylation_data[[1]])

  batch_id <- 1
  for(batch_id in 1:length(sample_sheet))
  {
    message("INFO: ", Sys.time(), " Working on batch:",batch_id)
    sample_sheet_local <- sample_sheet[[batch_id]]
    methylation_data_local <- stats::na.omit(methylation_data[[batch_id]])
    methylation_data_local <- methylation_data_local[rownames(methylation_data_local)%in%probes_to_preserve,]
    sample_sheet_local <- analyze_batch(methylation_data_local, sample_sheet_local, sliding_window_size, bonferroni_threshold,iqrTimes, batch_id)
    if(exists("sample_sheet_result"))
      sample_sheet_result <- plyr::rbind.fill(sample_sheet_result, sample_sheet_local)
    else
      sample_sheet_result <- sample_sheet_local
    utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  }

  sample_sheet <- sample_sheet_result
  utils::write.csv2(sample_sheet, file.path(ssEnv$result_folderData , "sample_sheet_result.csv"), row.names = F)
  message("INFO: ", Sys.time(), " Saving Sample Sheet with Results! ", Sys.time())


  if(length(sample_sheet$Sample_Group=="Reference")>0)
    sample_groups <- c("Reference","Control","Case")
  else
    sample_groups <- c("Control","Case")

  figures <- ssEnv$keys_figures[,1]
  markers <- ssEnv$keys_markers[,1]

  if(sum(ssEnv$keys_metaareas[,1]=="PROBE")==1)
  {
    subGroups <- c("")
    probes_prefix = "PROBES"
    mainGroupLabel =  "PROBE"
    subGroupLabel= "AREA"

    create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
    chrBed <- annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( inputBedDataFrame =  chrBed,markers = markers, file_prefix = "PROBES", groupColumnLabels = c("PROBE"))
  }

  if(sum(ssEnv$keys_metaareas[,1]=="CHR")==1)
  {
    subGroups <- c("CHR")
    probes_prefix = "PROBES_CHR_"
    mainGroupLabel =  "CHR"
    subGroupLabel= "AREA"

    create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
    chrBed <- annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( inputBedDataFrame =  chrBed,markers = markers, file_prefix = "CHR", groupColumnLabels = c("CHR"))
  }

  if(sum(ssEnv$keys_metaareas[,1]=="GENE")==1)
  {
    subGroups <- ssEnv$gene_subareas[,1]
    probes_prefix = "PROBES_Gene_"
    mainGroupLabel =  "GENE"
    subGroupLabel="AREA"

    create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

    geneBed <- annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( inputBedDataFrame =  geneBed,markers = markers, file_prefix = "GENE_AREA", groupColumnLabels = c("AREA"))
    create_heatmap( inputBedDataFrame =  geneBed,markers = markers, file_prefix = "GENE", groupColumnLabels = c("GENE"))
    create_heatmap( inputBedDataFrame =  geneBed,markers = markers, file_prefix = "GENE_PARTS", groupColumnLabels = c("GENE","AREA"))
  }


  if(sum(ssEnv$keys_metaareas[,1]=="ISLAND")==1)
  {
    probes_prefix <- "PROBES_Island_"
    subGroups <- ssEnv$island_subareas[,1]
    mainGroupLabel <- "ISLAND"
    subGroupLabel <- "RELATION_TO_CPGISLAND"
    create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

    islandBed <- annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( inputBedDataFrame =  islandBed,markers = markers, file_prefix = "RELATION_TO_CPGISLAND", groupColumnLabels = "RELATION_TO_CPGISLAND")
    create_heatmap( inputBedDataFrame =  islandBed,markers = markers, file_prefix = "ISLAND", groupColumnLabels = "ISLAND")
    create_heatmap( inputBedDataFrame =  islandBed,markers = markers, file_prefix = "ISLAND_PARTS", groupColumnLabels = c("ISLAND","RELATION_TO_CPGISLAND"))
  }

  if(sum(ssEnv$keys_metaareas[,1]=="DMR")==1)
  {
    subGroups <- c("DMR")
    probes_prefix = "PROBES_DMR_"
    mainGroupLabel =  "DMR"
    subGroupLabel="AREA"
    create_excel_pivot ( sample_groups =  sample_groups, figures =  figures,markers =  markers, subGroups =  subGroups, probes_prefix =   probes_prefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)

    dmrBed <- annotate_bed(sample_groups ,figures ,markers ,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)
    create_heatmap( inputBedDataFrame =  dmrBed,markers = markers, file_prefix = mainGroupLabel, groupColumnLabels = "DMR" )
  }

  # if (!is.null(geneBed))
  # {
  #    geneBed <- geneBed[,c("MAINGROUP","SAMPLEID","SUBAREA","VALUE","FIGURE","MARKER","SAMPLE_GROUP")]
  # }
  #
  # if (!is.null(dmrBed))
  # {
  #   dmrBed <- dmrBed[,c("MAINGROUP","SAMPLEID","SUBAREA","VALUE","FIGURE","MARKER","SAMPLE_GROUP")]
  # }
  #
  # if (!is.null(islandBed))
  # {
  #   islandBed <- islandBed[,c("MAINGROUP","SAMPLEID","SUBAREA","VALUE","FIGURE","MARKER","SAMPLE_GROUP")]
  # }
  #
  # totalBed <- rbind(geneBed, dmrBed, islandBed)
  # if (!is.null(totalBed) && nrow(totalBed)>0)
  #   create_heatmap( inputBedDataFrame =  totalBed,markers = markers, file_prefix = "GENOMIC_AREA", groupColumnLabels = 3)
  #
  # rm(populationControlRangeBetaValues)

  # message("Starting inference Analysis.")
  # inferenceAnalysis(ssEnv$result_folderData = ssEnv$result_folderData, ssEnv$session_folder= ssEnv$session_folder, inferenceDetails)
  # future::autoStopCluster(computationCluster)
  # doFuture::stopImplicitCluster()

  # geneontology_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData, fileName = fileName)
  # euristic_analysis_webgestalt(ssEnv$result_folderData = ssEnv$result_folderData)

  # if(length(methylation_data)>1)
  #   batch_correlation_check(ssEnv)

  close_env()
}
