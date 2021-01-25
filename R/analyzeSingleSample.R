#' Title
#'
#' @param values values of methylation
#' @param slidingWindowSize size of window sliding to calculate hypergeometric
#' @param resultFolder folder to store result files
#' @param thresholds threshold to use for comparison
#' @param comparison function to use for compare
#' @param sampleName name of sample to analyze
#' @param probeFeatures probes features probe, chr, start,end
#' @param subFileExtension extension to pre pend to file name
#' @param bonferroniThreshold bonferroni threshol to validate pVale
#'
#' @return list of lesion count  and probes count
#' @export
#'
analyzeSingleSample <- function(values, slidingWindowSize, resultFolder, thresholds, comparison, sampleName, subFileExtension, bonferroniThreshold = 0.05, probeFeatures) {

  #
  MUTATIONS <- NULL
  CHR <- NULL

  start_time_single_sample <- Sys.time()
  message(sampleName, " ", "... Sample analysis warmingUP ", Sys.time())
  result <- ""
  colnames(values) <- "VALUE"

  message(sampleName, " ", "Sample analysis WarmedUP ...", Sys.time())
  message(sampleName, " ", "Start sample analyze ", Sys.time())

  ### get probesOverThreshold ################################################################################################

  mutation <- as.numeric(comparison(values, thresholds))
  message(sampleName, " ", "Got probesOverThreshold ", Sys.time())
  ### get mutationAnnotatedSorted ################################################################################################
  if (!test_match_order(row.names(mutation), probeFeatures$PROBE)) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  mutationAnnotated <- data.frame(as.data.frame(probeFeatures), "MUTATIONS" = mutation, row.names = probeFeatures$PROBE)

  if (!test_match_order(row.names(mutationAnnotated), probeFeatures$PROBE)) {
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  }

  mutationAnnotatedSorted <- sortByCHRandSTART(mutationAnnotated)

  if (!test_match_order(mutationAnnotatedSorted$PROBE, row.names(mutationAnnotatedSorted))) {
    stop("Mutation annotation sorted is not coherent with probe informations order!")
  }

  result <- c("mutationCount" = sum(mutationAnnotatedSorted$MUTATIONS), "lesionCount" = 0, "probesCount" = 0)

  mutationAnnotatedSortedToSave <- subset(mutationAnnotatedSorted, MUTATIONS == 1)[, c("CHR", "START", "END")]
  message(sampleName, " ", "Got mutationAnnotatedSorted ", Sys.time())

  dumpSampleAsBedFile(
    dataToDump = mutationAnnotatedSortedToSave,
    fileExtension = paste(".", subFileExtension, ".MUTATIONS.bed", sep = ""),
    resultFolder = resultFolder,
    resultSubFolder = paste("MUTATIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  ### get lesion #################################################################################################
  lesionWeighted <- getLesions(mutationAnnotatedSorted = mutationAnnotatedSorted, slidingWindowSize = slidingWindowSize, sampleName = sampleName, bonferroniThreshold = bonferroniThreshold, probeFeatures = probeFeatures)

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  PROBES_Gene_5UTR <- get("PROBES_Gene_5UTR")
  gene_probes <- PROBES_Gene_5UTR
  gene_probes$CHR <- as.factor(gene_probes$CHR)
  gene_probes$END <- gene_probes$START
  gene_probes <- gene_probes[, c("PROBE","CHR","START","END","GENE")]
  gene_probes <- unique(gene_probes)

  mutationAnnotatedSortedGENE <- dplyr::inner_join(mutationAnnotatedSorted,gene_probes, by = c("CHR", "START","END","PROBE"))
  lesionWeighted <- getLesionsNew(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "GENE", mutationAnnotatedSorted = mutationAnnotatedSortedGENE)

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".GENE_5UTR.LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension,"GENE_5UTR", sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  PROBES_Gene_3UTR <- get("PROBES_Gene_3UTR")
  gene_probes <- PROBES_Gene_3UTR
  gene_probes$CHR <- as.factor(gene_probes$CHR)
  gene_probes$END <- gene_probes$START
  gene_probes <- gene_probes[, c("PROBE","CHR","START","END","GENE")]
  gene_probes <- unique(gene_probes)

  mutationAnnotatedSortedGENE <- dplyr::inner_join(mutationAnnotatedSorted,gene_probes, by = c("CHR", "START","END","PROBE"))
  lesionWeighted <- getLesionsNew(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "GENE", mutationAnnotatedSorted = mutationAnnotatedSortedGENE)

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".GENE_3UTR.LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension,"GENE_3UTR", sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  gene_probes <- get("PROBES_Gene_Body")
  gene_probes$CHR <- as.factor(gene_probes$CHR)
  gene_probes$END <- gene_probes$START
  gene_probes <- gene_probes[, c("PROBE","CHR","START","END","GENE")]
  gene_probes <- unique(gene_probes)

  mutationAnnotatedSortedGENE <- dplyr::inner_join(mutationAnnotatedSorted,gene_probes, by = c("CHR", "START","END","PROBE"))
  lesionWeighted <- getLesionsNew(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "GENE", mutationAnnotatedSorted = mutationAnnotatedSortedGENE)

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".GENE_BODY.LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension,"GENE_BODY", sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )


  dmr_probes <- get("PROBES_DMR_DMR")
  dmr_probes$CHR <- as.factor(dmr_probes$CHR)
  dmr_probes <- unique(dmr_probes)
  mutationAnnotatedSortedDMR <- dplyr::inner_join(mutationAnnotatedSorted,dmr_probes, by = c("CHR", "START","END","PROBE"))
  lesionWeighted <- getLesionsNew(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "DMR", mutationAnnotatedSorted = mutationAnnotatedSortedDMR)

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".DMR.LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension,"DMR_DMR", sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )


  lesionWeighted <- getLesionsNew(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "CHR", mutationAnnotatedSorted = mutationAnnotatedSorted)
  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".CHR.LESIONS.bed"),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension,"CHR", sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLENAME")
  )

  result["lesionCount"] <- dim(lesionWeighted)[1]
  result["probesCount"] <- dim(probeFeatures)[1]
  # if (result["lesionCount"] > dim(mutationAnnotatedSortedToSave)[1])
  # {
  #   ##
  #   lesionWeighted <- getLesions(mutationAnnotatedSorted = mutationAnnotatedSorted, slidingWindowSize = slidingWindowSize , sampleName = sampleName, probeFeatures =  probeFeatures)
  #   result["lesionCount"] <- dim(lesionWeighted)[1]
  # }


  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sampleName, " ", "Completed sample ", time_taken)
  return(result)
  # rm(list = ls())
}
