#' deltaSingleSample
#'
#' @param values values of methylation
#' @param highThresholds highest threshold to use for comparison
#' @param lowThresholds lowest threshold to use for comparison
#' @param sampleDetail details of sample to analyze
#' @param betaMedians median to use for calculation
#' @param probeFeatures genomic position of probes
#'
#' @return none

#'
deltaSingleSample <- function( values, highThresholds, lowThresholds, sampleDetail, betaMedians, probeFeatures) {

 message(sampleDetail$Sample_ID, " ", "DeltaSingleSample Sample analysis warmingUP ", Sys.time())

 # values <- data.frame("VALUE"=values)

 # browser()

 ### get probeFeatures ################################################################################################

 message(sampleDetail$Sample_ID, " ", "DeltaSingleSample Sample analysis WarmedUP ...", Sys.time())
 message(sampleDetail$Sample_ID, " ", "DeltaSingleSample Start sample analyze ", Sys.time())

 ### get probesOverThreshold ################################################################################################
 if (!test_match_order(row.names(values), probeFeatures$PROBE)) {
         stop("Wrong order matching Probes and Values!", Sys.time())
 }

 if (!test_match_order(row.names(highThresholds), probeFeatures$PROBE)) {
         stop("Wrong order matching Probes and highThresholds!", Sys.time())
 }

 if (!test_match_order(row.names(lowThresholds), probeFeatures$PROBE)) {
         stop("Wrong order matching Probes and lowThresholds!", Sys.time())
 }

 mutationAbove <- values > highThresholds
 mutationBelow <- values < lowThresholds
 mutation <- data.frame("MUTATIONS"= (mutationBelow & mutationAbove), row.names = probeFeatures$PROBE)
 colnames(mutation) <- "MUTATIONS"

 message(sampleDetail$Sample_ID, " ", "Got outliers ", Sys.time())

 ### get deltas #########################################################

 deltas <- data.frame("DELTA"= round(values - betaMedians,3))
 colnames(deltas) <- "DELTA"

 if (!test_match_order(row.names(mutation), probeFeatures$PROBE)) {
 stop("Wrong order matching Probes and Mutation!", Sys.time())
 }
 if (!test_match_order(row.names(deltas), probeFeatures$PROBE)) {
 stop("Wrong order matching Probes and Mutation!", Sys.time())
 }

 deltasAnnotated <- data.frame(as.data.frame(probeFeatures), deltas, "MUTATIONS" = mutation)

 deltasAnnotatedSorted <- sortByCHRandSTART(deltasAnnotated)

 deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, deltasAnnotatedSorted$MUTATIONS == 1)[, c("CHR", "START", "END", "DELTA")]

 folder2Save <- dir_check_and_create(resultFolderData,c(as.character(sampleDetail$Sample_Group),"DELTAS_METHYLATION"))
 # browser()
 result <- ""
 result <- result[-1]
 result["DELTA_AVG"] <- round(mean(deltas$DELTA),3)
 result["DELTA_VAR"] <- round(stats::var(deltas$DELTA),3)
 result["DELTA_MEDIAN"] <- round(stats::median(deltas$DELTA),3)

 dumpSampleAsBedFile( dataToDump = deltasAnnotatedSorted, fileName = file_path_build(folder2Save,c(as.character(sampleDetail$Sample_ID),"DELTAS","METHYLATION"),"bedgraph"))

 return(result)
}
