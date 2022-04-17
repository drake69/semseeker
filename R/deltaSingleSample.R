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

 message(sample_detail$Sample_ID, " ", "DeltaSingleSample Sample analysis warmingUP ", Sys.time())
 ### get probe_features ################################################################################################

 message(sample_detail$Sample_ID, " ", "DeltaSingleSample Sample analysis WarmedUP ...", Sys.time())
 message(sample_detail$Sample_ID, " ", "DeltaSingleSample Start sample analyze ", Sys.time())

 ### get probesOverThreshold ################################################################################################
 if (!test_match_order(row.names(values), probe_features$PROBE)) {
         stop("Wrong order matching Probes and Values!", Sys.time())
 }

 if (!test_match_order(row.names(high_thresholds), probe_features$PROBE)) {
         stop("Wrong order matching Probes and high_thresholds!", Sys.time())
 }

 if (!test_match_order(row.names(low_thresholds), probe_features$PROBE)) {
         stop("Wrong order matching Probes and low_thresholds!", Sys.time())
 }

 mutation_above <- values > high_thresholds
 mutation_below <- values < low_thresholds
 mutation <- data.frame("MUTATIONS"= (mutation_below & mutation_above), row.names = probe_features$PROBE) # nolint
 colnames(mutation) <- "MUTATIONS"

 message(sample_detail$Sample_ID, " ", "Got outliers ", Sys.time())

 ### get deltas #########################################################

 deltas <- data.frame("DELTA"= round(values - beta_medians,3))
 colnames(deltas) <- "DELTA"

 if (!test_match_order(row.names(mutation), probe_features$PROBE)) {
 stop("Wrong order matching Probes and Mutation!", Sys.time())
 }
 if (!test_match_order(row.names(deltas), probe_features$PROBE)) {
 stop("Wrong order matching Probes and Mutation!", Sys.time())
 }

 deltasAnnotated <- data.frame(as.data.frame(probe_features), deltas, "MUTATIONS" = mutation)
 deltasAnnotatedSorted <- sort_by_chr_and_start(deltasAnnotated)
 deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, deltasAnnotatedSorted$MUTATIONS == 1)[, c("CHR", "START", "END", "DELTA")]

 result <- ""
 result <- result[-1]
 result["DELTA_AVG"] <- round(mean(deltas$DELTA),3)
 result["DELTA_VAR"] <- round(stats::var(deltas$DELTA),3)
 result["DELTA_MEDIAN"] <- round(stats::median(deltas$DELTA),3)

 message("############# SEARCH")
 message("############# SEARCH",search())
 message("############# LS",ls())
 # message("############# envir$result_folderData:", envir$result_folderData)

 folder_to_save <- dir_check_and_create(envir$result_folderData,c(as.character(sample_detail$Sample_Group),"DELTAS_METHYLATION"))
 dump_sample_as_bed_file(dataToDump = deltasAnnotatedSorted, fileName = file_path_build(folder_to_save,c(as.character(sample_detail$Sample_ID),"DELTAS","METHYLATION"),"bedgraph"))
 return(result)
}
