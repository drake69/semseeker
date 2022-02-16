#' analyzeSingleSample
#'
#' @param values values of methylation
#' @param slidingWindowSize size of window sliding to calculate hypergeometric
#' @param thresholds threshold to use for comparison
#' @param sampleDetail details of the sample to analyze
#' @param figure which figure's of sasmple will be analized HYPO or HYPER
#' @param bonferroniThreshold bonferroni threshold to validate pValue
#' @param probeFeatures probes details to be used
#' @param envir environment to get globals
#' @return list of lesion count  and probes count
#'
#'
analyzeSingleSample <- function(envir, values, slidingWindowSize, thresholds, figure, sampleDetail, bonferroniThreshold = 0.05, probeFeatures) {

  # browser()
  start_time_single_sample <- Sys.time()
  message(sampleDetail$Sample_ID, " ", "SingleSample Sample analysis warmingUP ", Sys.time())
  result <- ""
  result <- result[-1]

  # colnames(values) <- "VALUE"

  message(sampleDetail$Sample_ID, " ", "SingleSample Sample analysis WarmedUP ...", Sys.time())
  message(sampleDetail$Sample_ID, " ", "SingleSample Start sample analyze ", Sys.time())

  mutationAnnotatedSorted <- getMutations(values, figure,thresholds, probeFeatures, sampleDetail$Sample_ID)
  mutationAnnotatedSortedToSave <- subset(mutationAnnotatedSorted, mutationAnnotatedSorted$MUTATIONS == 1)[, c("CHR", "START", "END")]

  message("############# SEARCH")
  message("############# SEARCH",search())
  message("############# LS",ls())
  # browser()
  message("############# envir$resultFolderData:", envir$resultFolderData)
  folder2Save <- dir_check_and_create(envir$resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  dumpSampleAsBedFile(
    dataToDump = mutationAnnotatedSortedToSave,
    fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"MUTATIONS",figure),"bed")
  )
  result[paste("MUTATIONS_", figure, sep="")] <- if (!is.null(mutationAnnotatedSortedToSave)) dim(mutationAnnotatedSortedToSave)[1] else 0

  lesionWeighted <- getLesions(bonferroniThreshold = bonferroniThreshold, slidingWindowSize = slidingWindowSize, grouping_column = "CHR", mutationAnnotatedSorted = mutationAnnotatedSorted)
  folder2Save <- dir_check_and_create(envir$resultFolderData,c(as.character(sampleDetail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileName = file_path_build(folder2Save,c(sampleDetail$Sample_ID,"LESIONS",figure),"bed")
  )

  result[paste("LESIONS_", figure, sep="")] <- if (!is.null(lesionWeighted)) dim(lesionWeighted)[1] else 0
  if(figure=="HYPER")
     result["PROBES_COUNT"] <- dim(probeFeatures)[1]

  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sampleDetail$Sample_ID, " ", "Completed sample ", time_taken)
  return(result)
}



