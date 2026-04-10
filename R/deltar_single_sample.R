#' deltar_single_sample
#'
#' @param values data.frame of methylation values with columns CHR, START, END and signal value in column 4
#' @param thresholds data.frame of signal thresholds (from signal_range_values) with columns CHR, START, END,
#'   signal_superior_thresholds, signal_inferior_thresholds, signal_median_values
#' @param sample_detail named list/row with at least Sample_ID and Sample_Group fields
#' @return invisibly NULL; relative-delta results are written as bedgraph.gz files under the session data folder
#'
deltar_single_sample <- function ( values, thresholds, sample_detail) {

  ssEnv <- get_session_info()
  result <- data.frame()
  result <- result[-1]

  values <- sort_by_chr_and_start(values)
  thresholds <- sort_by_chr_and_start(thresholds)

  high_thresholds <-  thresholds$signal_superior_thresholds
  signal_medians <- thresholds$signal_median_values
  low_thresholds <- thresholds$signal_inferior_thresholds

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
  deltar_hyper <- data.frame("CHR"= thresholds$CHR,"START"= thresholds$START, "END"=thresholds$END, "DELTA"= (values[,4] - high_thresholds)/dividend)
  deltar_hyper <- sort_by_chr_and_start(deltar_hyper)
  # save only deltar of epimutation
  deltar_hyper <- subset(deltar_hyper, deltar_hyper$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltar_hyper, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPER"),"bedgraph", add_gz=TRUE))

  ### get deltar_hypo HYPER #########################################################
  # if(length(low_thresholds)!=length(high_thresholds) || nrow(values[,4])!=length(high_thresholds))
  # {
  #   stop("ERROR: I'm stopping here the values and thresholds have different number of rows!")
  # }
  deltar_hypo <- data.frame("CHR"= thresholds$CHR,"START"= thresholds$START, "END"=thresholds$END, "DELTA"=  (low_thresholds - values[,4])/dividend)
  deltar_hypo <- sort_by_chr_and_start(deltar_hypo)
  # save only deltar of epimutation
  deltar_hypo <- subset(deltar_hypo, deltar_hypo$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltar_hypo, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPO"),"bedgraph", add_gz=TRUE))


  # if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAR_BOTH"))
  # {
  #   ### get deltar BOTH #########################################################
  #   deltarAnnotated_both <- rbind(deltar_hypo, deltar_hyper)
  #   deltarAnnotated_bothSorted <- sort_by_chr_and_start(deltarAnnotated_both)
  #   deltarAnnotated_bothSorted <- subset(deltarAnnotated_bothSorted, deltarAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]
  #
  #   folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_BOTH"))
  #   dump_sample_as_bed_file(data_to_dump = deltarAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","BOTH"),"bedgraph", add_gz=TRUE))
  # }
  #
  # if (any(ssEnv$keys_markers_figures$COMBINED=="DELTAR_BOTHSUM"))
  # {
  #   ### get deltar BOTHSUM #########################################################
  #   deltarAnnotated_hypoSorted_sign <- deltar_hypo
  #   deltarAnnotated_hypoSorted_sign$DELTA <- (-1 * deltarAnnotated_hypoSorted_sign$DELTA)
  #   deltarAnnotated_both_sum <- rbind(deltar_hyper, deltarAnnotated_hypoSorted_sign)
  #   deltarAnnotated_both_sumSorted <- sort_by_chr_and_start(deltarAnnotated_both_sum)
  #   deltarAnnotated_both_sumSorted <- subset(deltarAnnotated_both_sumSorted, deltarAnnotated_both_sumSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]
  #
  #   folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_BOTHSUM"))
  #   dump_sample_as_bed_file(data_to_dump = deltarAnnotated_both_sumSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","BOTHSUM"),"bedgraph", add_gz=TRUE))
  # }

  ### get deltar from medians #########################################################

  # deltar <- data.frame("DELTA"= values[,4] - signal_medians)
  # colnames(deltar) <- "DELTA"

  deltar_to_check <- c(deltar_hypo$DELTA, deltar_hyper$DELTA)
  if (length(deltar_to_check)>0)
    if( (min(deltar_to_check)<0))
    {
      log_event(min(deltar_hypo$DELTA))
      log_event(min(deltar_hyper$DELTA))
      stop("ERROR: I'm stopping here the deltar have negative values!")
    }
}





