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
delta_single_sample <- function ( values, thresholds, sample_detail) {

  ssEnv <- get_session_info()

  values <- sort_by_chr_and_start(values)
  thresholds <- sort_by_chr_and_start(thresholds)

  high_thresholds <-  thresholds$signal_superior_thresholds
  signal_medians <- thresholds$signal_median_values
  low_thresholds <- thresholds$signal_inferior_thresholds

  if(any(high_thresholds<low_thresholds))
    stop("ERROR: I'm stopping here the some high have values less than low!")

  ### get deltas HYPER #########################################################
  deltas_hyper <- data.frame("CHR"= thresholds$CHR,"START"= thresholds$START, "END"=thresholds$END, "DELTA"= values[,4] - high_thresholds)
  deltas_hyper <- sort_by_chr_and_start(deltas_hyper)
  # save only deltas of epimutation
  deltas_hyper <- subset(deltas_hyper, deltas_hyper$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltas_hyper, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPER"),"bedgraph", add_gz=TRUE))

  ### get deltas HYPO #########################################################
  deltas_hypo <- data.frame("CHR"= thresholds$CHR,"START"= thresholds$START, "END"=thresholds$END, "DELTA"=  low_thresholds - values[,4])
  deltas_hypo <- sort_by_chr_and_start(deltas_hypo)
  # save only deltas of epimutation
  deltas_hypo <- subset(deltas_hypo, deltas_hypo$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltas_hypo, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPO"),"bedgraph", add_gz=TRUE))

  # if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAS_BOTH"))
  # {
  #   ### get deltas BOTH #########################################################
  #   deltasAnnotated_both <- rbind(deltasAnnotated_hypoSorted, deltas_hyper)
  #   deltasAnnotated_bothSorted <- sort_by_chr_and_start(deltasAnnotated_both)
  #   deltasAnnotated_bothSorted <- subset(deltasAnnotated_bothSorted, deltasAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]
  #
  #   folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTH"))
  #   dump_sample_as_bed_file(data_to_dump = deltasAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTH"),"bedgraph", add_gz=TRUE))
  # }
  # ### get deltas BOTHSUM #########################################################
  #
  # if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAS_BOTHSUM"))
  # {
  #   deltasAnnotated_hypoSorted_sign <- deltasAnnotated_hypoSorted
  #   deltasAnnotated_hypoSorted_sign$DELTA <- (-1 * deltasAnnotated_hypoSorted_sign$DELTA)
  #   deltasAnnotated_both_sum <- rbind(deltas_hyper, deltasAnnotated_hypoSorted_sign)
  #   deltasAnnotated_both_sumSorted <- sort_by_chr_and_start(deltasAnnotated_both_sum)
  #   deltasAnnotated_both_sumSorted <- subset(deltasAnnotated_both_sumSorted, deltasAnnotated_both_sumSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]
  #
  #   folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTHSUM"))
  #   dump_sample_as_bed_file(data_to_dump = deltasAnnotated_both_sumSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTHSUM"),"bedgraph", add_gz=TRUE))
  # }

  ### get deltas from medians #########################################################

  # deltas <- data.frame("DELTA"= values - signal_medians, row.names = probe_features$PROBE)
  # colnames(deltas) <- "DELTA"

  deltas_to_check <- c(deltas_hypo$DELTA, deltas_hyper$DELTA)
  if (length(deltas_to_check) > 0)
    if( (min(deltas_to_check)<0))
    {
      log_event(min(deltasAnnotated_hypoSorted$DELTA))
      log_event(min(deltas_hyper$DELTA))
      stop("ERROR: I'm stopping here the deltas have negative values!")
    }
}






