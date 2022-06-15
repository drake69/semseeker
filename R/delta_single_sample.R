#' delta_single_sample
#'
#' @param values values of methylation
#' @param high_thresholds highest threshold to use for comparison
#' @param low_thresholds lowest threshold to use for comparison
#' @param sample_detail details of sample to analyze
#' @param beta_medians median to use for calculation
#' @param probe_features genomic position of probes
#' @param envir environment to get globals
#' @return summary detail about the analysis
#'
delta_single_sample <- function (envir, values, high_thresholds, low_thresholds, sample_detail, beta_medians, probe_features) {

 if (!test_match_order(row.names(values), probe_features$PROBE)) {
         stop("Wrong order matching Probes and Values!", Sys.time())
 }

 if (!test_match_order(row.names(high_thresholds), probe_features$PROBE)) {
         stop("Wrong order matching Probes and high_thresholds!", Sys.time())
 }

 if (!test_match_order(row.names(low_thresholds), probe_features$PROBE)) {
         stop("Wrong order matching Probes and low_thresholds!", Sys.time())
 }

 ### get deltas HYPER #########################################################
 deltas <- data.frame("DELTA"= round(values - high_thresholds,5))
 colnames(deltas) <- "DELTA"

 if (!test_match_order(row.names(deltas), probe_features$PROBE)) {
 stop("Wrong order matching Probes and Mutation!", Sys.time())
 }

 deltasAnnotated_hyper <- data.frame(as.data.frame(probe_features), deltas)
 deltasAnnotated_hyperSorted <- sort_by_chr_and_start(deltasAnnotated_hyper)
 deltasAnnotated_hyperSorted <- subset(deltasAnnotated_hyperSorted, deltasAnnotated_hyperSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

 folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPER"))
 dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hyperSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPER"),"bedgraph"))

 ### get deltas HYPER #########################################################
 deltas <- data.frame("DELTA"= round( - values + low_thresholds,5))
 colnames(deltas) <- "DELTA"

 if (!test_match_order(row.names(deltas), probe_features$PROBE)) {
   stop("Wrong order matching Probes and Mutation!", Sys.time())
 }

 deltasAnnotated_hypo <- data.frame(as.data.frame(probe_features), deltas)
 deltasAnnotated_hypoSorted <- sort_by_chr_and_start(deltasAnnotated_hypo)
 deltasAnnotated_hypoSorted <- subset(deltasAnnotated_hypoSorted, deltasAnnotated_hypoSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

 folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_HYPO"))
 dump_sample_as_bed_file(data_to_dump = deltasAnnotated_hypoSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","HYPO"),"bedgraph"))


 ### get deltas BOTH #########################################################
 deltasAnnotated_both <- rbind(deltasAnnotated_hyper, deltasAnnotated_hypo)
 deltasAnnotated_bothSorted <- sort_by_chr_and_start(deltasAnnotated_both)
 deltasAnnotated_bothSorted <- subset(deltasAnnotated_bothSorted, deltasAnnotated_bothSorted$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

 folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_BOTH"))
 dump_sample_as_bed_file(data_to_dump = deltasAnnotated_bothSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","BOTH"),"bedgraph"))


 ### get deltas from medians #########################################################

 deltas <- data.frame("DELTA"= round(values - beta_medians,5))
 colnames(deltas) <- "DELTA"

 if (!test_match_order(row.names(deltas), probe_features$PROBE)) {
   stop("Wrong order matching Probes and Mutation!", Sys.time())
 }

 deltasAnnotated <- data.frame(as.data.frame(probe_features), deltas)
 deltasAnnotatedSorted <- sort_by_chr_and_start(deltasAnnotated)
 deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, deltasAnnotatedSorted$DELTA != 0)[, c("CHR", "START", "END", "DELTA")]

 result <- ""
 result <- result[-1]
 result["DELTA_AVG"] <- round(mean(deltasAnnotatedSorted$DELTA),5)
 result["DELTA_VAR"] <- round(stats::var(deltasAnnotatedSorted$DELTA),5)
 result["DELTA_MEDIAN"] <- round(stats::median(deltasAnnotatedSorted$DELTA),5)


 return(result)
}






