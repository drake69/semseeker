#' delta_single_sample
#'
#' @param values values of methylation
#' @param high_thresholds highest threshold to use for comparison
#' @param low_thresholds lowest threshold to use for comparison
#' @param sample_detail details of sample to analyze
#' @param beta_medians median to use for calculation
#' @param probe_features genomic position of probe_features
#' @return summary detail about the analysis
#'
deltar_single_sample <- function ( values, high_thresholds, low_thresholds, sample_detail, beta_medians, probe_features) {

  ssEnv <- get_session_info()


  dividend <- high_thresholds - low_thresholds
  names(dividend) <- "DIVIDEND"
  dividend[dividend==0] <- 0.000000001

  if(any(dividend<0))
    stop("ERROR: I'm stopping here the dividend have negative values!")


  ### get deltar HYPER #########################################################
  deltar_hyper <- data.frame("DELTA"= (values - high_thresholds)/dividend, row.names = probe_features$PROBE)
  colnames(deltar_hyper) <- "DELTA"
  deltarAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltar_hyper)
  deltarAnnotated_hyperSorted <- sort_by_chr_and_start(deltarAnnotated_hyper)
  # save only deltar of epimutation
  deltarAnnotated_hyperSorted <- subset(deltarAnnotated_hyperSorted, deltarAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltarAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPER"),"bedgraph"))

  ### get deltar_hypo HYPER #########################################################
  deltar_hypo <- data.frame("DELTA"=  (low_thresholds - values)/dividend, row.names = probe_features$PROBE)
  colnames(deltar_hypo) <- "DELTA"
  deltarAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltar_hypo, row.names = probe_features$PROBE)
  deltarAnnotated_hypoSorted <- sort_by_chr_and_start(deltarAnnotated_hypo)
  # save only deltar of epimutation
  deltarAnnotated_hypoSorted <- subset(deltarAnnotated_hypoSorted, deltarAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltarAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","HYPO"),"bedgraph"))


  ### get deltar BOTH #########################################################
  deltarAnnotated_both <- rbind(deltarAnnotated_hypoSorted, deltarAnnotated_hyperSorted)
  deltarAnnotated_bothSorted <- sort_by_chr_and_start(deltarAnnotated_both)
  deltarAnnotated_bothSorted <- subset(deltarAnnotated_bothSorted, deltarAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAR_BOTH"))
  dump_sample_as_bed_file(data_to_dump = deltarAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAR","BOTH"),"bedgraph"))


  ### get deltar from medians #########################################################

  # deltar <- data.frame("DELTA"= values - beta_medians, row.names = probe_features$PROBE)
  # colnames(deltar) <- "DELTA"

  result <- ""
  result <- result[-1]
  if(mean(deltarAnnotated_hypoSorted$DELTA)<0 |
      mean(deltarAnnotated_hyperSorted$DELTA) <0 |
      mean(deltarAnnotated_bothSorted$DELTA) <0 )
  {
    message(mean(deltarAnnotated_hypoSorted$DELTA))
    message(mean(deltarAnnotated_hyperSorted$DELTA))
    message(mean(deltarAnnotated_bothSorted$DELTA))
    stop("ERROR: I'm stopping here the deltar have negative values!")
  }
  result["DELTAR_HYPO"] <- mean(deltarAnnotated_hypoSorted$DELTA)
  result["DELTAR_HYPER"] <- mean(deltarAnnotated_hyperSorted$DELTA)
  result["DELTAR_BOTH"] <- mean(deltarAnnotated_bothSorted$DELTA)


  return(result)
}






