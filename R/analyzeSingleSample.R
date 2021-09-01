#' analyzeSingleSample
#'
#' @param values values of methylation
#' @param slidingWindowSize size of window sliding to calculate hypergeometric
#' @param resultFolder folder to store result files
#' @param thresholds threshold to use for comparison
#' @param comparison function to use for compare
#' @param sampleName name of sample to analyze
#' @param probeFeatures probes features probe, chr, start,end
#' @param subFileExtension extension to pre pend to file name
#' @param bonferroniThreshold bonferroni threshold to validate pValue
#'
#' @return list of lesion count  and probes count
#'
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

  mutationAnnotatedSorted <- getMutations(values, comparison,thresholds, probeFeatures, sampleName)

  mutationAnnotatedSortedToSave <- subset(mutationAnnotatedSorted, MUTATIONS == 1)[, c("CHR", "START", "END")]
  dumpSampleAsBedFile(
    dataToDump = mutationAnnotatedSortedToSave,
    fileExtension = paste(".", subFileExtension, ".MUTATIONS.bed", sep = ""),
    resultFolder = resultFolder,
    resultSubFolder = paste("MUTATIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID")
  )
  result["mutationCount"] <- if (!is.null(mutationAnnotatedSortedToSave)) dim(mutationAnnotatedSortedToSave)[1] else 0


  lesionWeighted <- getLesions(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "CHR", mutationAnnotatedSorted = mutationAnnotatedSorted)
  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension, ".LESIONS.bed", sep =""),
    resultFolder = resultFolder,
    resultSubFolder = paste0("LESIONS","_", subFileExtension, sep = ""),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID")
  )

  result["lesionCount"] <- if (!is.null(lesionWeighted)) dim(lesionWeighted)[1] else 0
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
}
