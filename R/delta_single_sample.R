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
delta_single_sample <- function ( values, high_thresholds, low_thresholds, sample_detail, beta_medians, probe_features) {

  ssEnv <- get_session_info()

  if(any(high_thresholds<low_thresholds))
    stop("ERROR: I'm stopping here the some high have values less than low!")


  ### get deltas HYPER #########################################################
  deltas_hyper <- data.frame("DELTA"= values - high_thresholds, row.names = probe_features$PROBE)
  colnames(deltas_hyper) <- "DELTA"
  deltasAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltas_hyper)
  deltasAnnotated_hyperSorted <- sort_by_chr_and_start(deltasAnnotated_hyper)
  # save only deltas of epimutation
  deltasAnnotated_hyperSorted <- subset(deltasAnnotated_hyperSorted, deltasAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPER"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPER"),"bedgraph"))

  ### get deltas_hypo HYPER #########################################################
  deltas_hypo <- data.frame("DELTA"=  low_thresholds - values, row.names = probe_features$PROBE)
  colnames(deltas_hypo) <- "DELTA"
  deltasAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltas_hypo, row.names = probe_features$PROBE)
  deltasAnnotated_hypoSorted <- sort_by_chr_and_start(deltasAnnotated_hypo)
  # save only deltas of epimutation
  deltasAnnotated_hypoSorted <- subset(deltasAnnotated_hypoSorted, deltasAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPO"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPO"),"bedgraph"))


  ### get deltas BOTH #########################################################
  deltasAnnotated_both <- rbind(deltasAnnotated_hypoSorted, deltasAnnotated_hyperSorted)
  deltasAnnotated_bothSorted <- sort_by_chr_and_start(deltasAnnotated_both)
  deltasAnnotated_bothSorted <- subset(deltasAnnotated_bothSorted, deltasAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTH"))
  dump_sample_as_bed_file(data_to_dump = deltasAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTH"),"bedgraph"))


  ### get deltas from medians #########################################################

  # deltas <- data.frame("DELTA"= values - beta_medians, row.names = probe_features$PROBE)
  # colnames(deltas) <- "DELTA"

  result <- ""
  result <- result[-1]
  if(mean(deltasAnnotated_hypoSorted$DELTA)<0 |
      mean(deltasAnnotated_hyperSorted$DELTA) <0 |
      mean(deltasAnnotated_bothSorted$DELTA) <0 )
  {
    message(mean(deltasAnnotated_hypoSorted$DELTA))
    message(mean(deltasAnnotated_hyperSorted$DELTA))
    message(mean(deltasAnnotated_bothSorted$DELTA))
    stop("ERROR: I'm stopping here the deltas have negative values!")
  }
  result["DELTAS_HYPO"] <- mean(deltasAnnotated_hypoSorted$DELTA)
  result["DELTAS_HYPER"] <- mean(deltasAnnotated_hyperSorted$DELTA)
  result["DELTAS_BOTH"] <- mean(deltasAnnotated_bothSorted$DELTA)


  return(result)
}






