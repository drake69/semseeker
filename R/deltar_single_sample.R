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
deltar_single_sample <- function ( values, high_thresholds, low_thresholds, sample_detail, signal_medians, probe_features) {

  ssEnv <- get_session_info()
  result <- data.frame()
  result <- result[-1]

  chk1 <- any(is.na(high_thresholds))
  chk2 <- any(is.na(low_thresholds))
  if(chk1 | chk2)
  {

  }

  dividend <- high_thresholds - low_thresholds
  names(dividend) <- "DIVIDEND"

  chk3 <- any(is.na(dividend))
  chk4 <- any(dividend<0)
  if(chk3 | chk4)
  {

    stop("ERROR: I'm stopping here the dividend have negative values!")
  }

  chk0 <- any(dividend==0)
  if(chk0)
    dividend[dividend==0] <- 0.000000001

  ### get deltar HYPER #########################################################
  deltar_hyper <- data.frame("DELTA"= (values - high_thresholds)/dividend, row.names = probe_features$PROBE)
  colnames(deltar_hyper) <- "DELTA"
  deltarAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltar_hyper)
  deltarAnnotated_hyperSorted <- sort_by_chr_and_start(deltarAnnotated_hyper)
  # save only deltar of epimutation
  deltarAnnotated_hyperSorted <- subset(deltarAnnotated_hyperSorted, deltarAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltarAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPER"),"bedgraph", add_gz=TRUE))

  ### get deltar_hypo HYPER #########################################################
  deltar_hypo <- data.frame("DELTA"=  (low_thresholds - values)/dividend, row.names = probe_features$PROBE)
  colnames(deltar_hypo) <- "DELTA"
  deltarAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltar_hypo, row.names = probe_features$PROBE)
  deltarAnnotated_hypoSorted <- sort_by_chr_and_start(deltarAnnotated_hypo)
  # save only deltar of epimutation
  deltarAnnotated_hypoSorted <- subset(deltarAnnotated_hypoSorted, deltarAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltarAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPO"),"bedgraph", add_gz=TRUE))


  if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAR_BOTH"))
  {
    ### get deltar BOTH #########################################################
    deltarAnnotated_both <- rbind(deltarAnnotated_hypoSorted, deltarAnnotated_hyperSorted)
    deltarAnnotated_bothSorted <- sort_by_chr_and_start(deltarAnnotated_both)
    deltarAnnotated_bothSorted <- subset(deltarAnnotated_bothSorted, deltarAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_BOTH"))
    dump_sample_as_bed_file(data_to_dump = deltarAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","BOTH"),"bedgraph", add_gz=TRUE))
    result <- data.frame_add.column(result, "DELTAR_BOTH",mean(deltarAnnotated_bothSorted$DELTA, na.rm = TRUE))
  }

  if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAR_BOTHSUM"))
  {
    ### get deltar BOTHSUM #########################################################
    deltarAnnotated_hypoSorted_sign <- deltarAnnotated_hypoSorted
    deltarAnnotated_hypoSorted_sign$DELTA <- (-1 * deltarAnnotated_hypoSorted_sign$DELTA)
    deltarAnnotated_both_sum <- rbind(deltarAnnotated_hyperSorted, deltarAnnotated_hypoSorted_sign)
    deltarAnnotated_both_sumSorted <- sort_by_chr_and_start(deltarAnnotated_both_sum)
    deltarAnnotated_both_sumSorted <- subset(deltarAnnotated_both_sumSorted, deltarAnnotated_both_sumSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]

    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_BOTHSUM"))
    dump_sample_as_bed_file(data_to_dump = deltarAnnotated_both_sumSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","BOTHSUM"),"bedgraph", add_gz=TRUE))
    result <- data.frame_add.column(result, "DELTAR_BOTHSUM",mean(deltarAnnotated_both_sumSorted$DELTA, na.rm = TRUE))
  }

  ### get deltar from medians #########################################################

  # deltar <- data.frame("DELTA"= values - signal_medians, row.names = probe_features$PROBE)
  # colnames(deltar) <- "DELTA"

  deltar_to_check <- c(deltarAnnotated_hypoSorted$DELTA, deltarAnnotated_hyperSorted$DELTA)
  if (length(deltar_to_check)>0)
    if( (min(deltar_to_check)<0))
    {
      log_event(min(deltarAnnotated_hypoSorted$DELTA))
      log_event(min(deltarAnnotated_hyperSorted$DELTA))
      stop("ERROR: I'm stopping here the deltar have negative values!")
    }
  result <- data.frame_add.column(result, "DELTAR_HYPO",mean(deltarAnnotated_hypoSorted$DELTA, na.rm = TRUE))
  result <- data.frame_add.column(result, "DELTAR_HYPER",mean(deltarAnnotated_hyperSorted$DELTA, na.rm = TRUE))
  # fill NA with 0
  result[is.na(result)] <- 0

  return(result)
}





