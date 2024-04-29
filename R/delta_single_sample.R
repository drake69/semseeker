#' delta_single_sample
#'
#' @param values values of methylation
#' @param high_thresholds highest threshold to use for comparison
#' @param low_thresholds lowest threshold to use for comparison
#' @param sample_detail details of sample to analyze
#' @param signal_medians median to use for calculation
#' @param probe_features genomic position of probe_features
#' @return summary detail about the analysis
#'
delta_single_sample <- function ( values, high_thresholds, low_thresholds, sample_detail, signal_medians, probe_features) {

  ssEnv <- get_session_info()

  if(any(high_thresholds<low_thresholds))
    stop("ERROR: I'm stopping here the some high have values less than low!")

  result <- data.frame()
  result <- result[-1]

  ### get deltas HYPER #########################################################
  deltas_hyper <- data.frame("DELTA"= values - high_thresholds, row.names = probe_features$PROBE)
  colnames(deltas_hyper) <- "DELTA"
  deltasAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltas_hyper)
  deltasAnnotated_hyperSorted <- sort_by_chr_and_start(deltasAnnotated_hyper)
  # save only deltas of epimutation
  deltasAnnotated_hyperSorted <- subset(deltasAnnotated_hyperSorted, deltasAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPER"),"bedgraph", add_gz=TRUE))

  ### get deltas_hypo HYPER #########################################################
  deltas_hypo <- data.frame("DELTA"=  low_thresholds - values, row.names = probe_features$PROBE)
  colnames(deltas_hypo) <- "DELTA"
  deltasAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltas_hypo, row.names = probe_features$PROBE)
  deltasAnnotated_hypoSorted <- sort_by_chr_and_start(deltasAnnotated_hypo)
  # save only deltas of epimutation
  deltasAnnotated_hypoSorted <- subset(deltasAnnotated_hypoSorted, deltasAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPO"),"bedgraph", add_gz=TRUE))

  if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAS_BOTH"))
  {
    ### get deltas BOTH #########################################################
    deltasAnnotated_both <- rbind(deltasAnnotated_hypoSorted, deltasAnnotated_hyperSorted)
    deltasAnnotated_bothSorted <- sort_by_chr_and_start(deltasAnnotated_both)
    deltasAnnotated_bothSorted <- subset(deltasAnnotated_bothSorted, deltasAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTH"))
    dump_sample_as_bed_file(data_to_dump = deltasAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTH"),"bedgraph", add_gz=TRUE))
    result <- data.frame_add.column(result, "DELTAS_BOTH",mean(deltasAnnotated_bothSorted$DELTA, na.rm = TRUE))
  }
  ### get deltas BOTHSUM #########################################################

  if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAS_BOTHSUM"))
  {
    deltasAnnotated_hypoSorted_sign <- deltasAnnotated_hypoSorted
    deltasAnnotated_hypoSorted_sign$DELTA <- (-1 * deltasAnnotated_hypoSorted_sign$DELTA)
    deltasAnnotated_both_sum <- rbind(deltasAnnotated_hyperSorted, deltasAnnotated_hypoSorted_sign)
    deltasAnnotated_both_sumSorted <- sort_by_chr_and_start(deltasAnnotated_both_sum)
    deltasAnnotated_both_sumSorted <- subset(deltasAnnotated_both_sumSorted, deltasAnnotated_both_sumSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]

    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTHSUM"))
    dump_sample_as_bed_file(data_to_dump = deltasAnnotated_both_sumSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTHSUM"),"bedgraph", add_gz=TRUE))
    result <- data.frame_add.column(result, "DELTAS_BOTHSUM",mean(deltasAnnotated_both_sumSorted$DELTA, na.rm = TRUE))
  }

  ### get deltas from medians #########################################################

  # deltas <- data.frame("DELTA"= values - signal_medians, row.names = probe_features$PROBE)
  # colnames(deltas) <- "DELTA"

  deltas_to_check <- c(deltasAnnotated_hypoSorted$DELTA, deltasAnnotated_hyperSorted$DELTA)
  if (length(deltas_to_check) > 0)
    if( (min(deltas_to_check)<0))
    {
      log_event(min(deltasAnnotated_hypoSorted$DELTA))
      log_event(min(deltasAnnotated_hyperSorted$DELTA))
      stop("ERROR: I'm stopping here the deltas have negative values!")
    }

  result <- data.frame_add.column(result, "DELTAS_HYPO",mean(deltasAnnotated_hypoSorted$DELTA, na.rm = TRUE))
  result <- data.frame_add.column(result, "DELTAS_HYPER",mean(deltasAnnotated_hyperSorted$DELTA, na.rm = TRUE))
  # fill NA with 0
  result[is.na(result)] <- 0

  return(result)
}






