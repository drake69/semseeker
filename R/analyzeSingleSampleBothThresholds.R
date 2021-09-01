#' analyzeSingleSampleBothThresholds
#'
#' @param values values of methylation
#' @param slidingWindowSize size of window sliding to calculate hypergeometric
#' @param resultFolder folder to store result files
#' @param betaSuperiorThresholds threshold for hyper methylation
#' @param betaInferiorThresholds threshold for hypo methylation
#' @param sampleName name of sample to analyze
#' @param probeFeatures probes features probe, chr, start,end
#' @param subFileExtension extension to pre pend to file name
#'
#' @return list of lesion count  and probes count
#'
#'
analyzeSingleSampleBothThresholds <- function(values, resultFolder, betaSuperiorThresholds, betaInferiorThresholds, sampleName, probeFeatures) {

  #
  MUTATIONS <- NULL
  CHR <- NULL

  subFileExtension <- "BOTH"
  start_time_single_sample <- Sys.time()
  message(sampleName, " ", "... Sample analysis warmingUP ", Sys.time())
  colnames(values) <- "VALUE"

  mutationAnnotatedSortedHyper <- getMutations(values =  values,comparison =  `>`,thresholds =  betaSuperiorThresholds,probeFeatures =  probeFeatures,sampleName =  sampleName)
  mutationAnnotatedSortedHypo <- getMutations(values =  values,comparison =  `<`,thresholds =  betaInferiorThresholds,probeFeatures =  probeFeatures,sampleName =  sampleName)
  mutationAnnotatedSorted <- rbind(mutationAnnotatedSortedHypo, mutationAnnotatedSortedHyper)
  mutationAnnotatedSorted <- sortByCHRandSTART(mutationAnnotatedSorted)
  mutationAnnotatedSortedToSave <- subset(mutationAnnotatedSorted, MUTATIONS == 1)[, c("CHR", "START", "END")]


  dumpSampleAsBedFile(
    dataToDump = mutationAnnotatedSortedToSave,
    fileExtension = paste(".", subFileExtension, ".MUTATIONS.bed", sep = ""),
    resultFolder = resultFolder,
    resultSubFolder = paste("MUTATIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID")
  )

  if (file.exists(paste(resultFolder, "/LESIONS_HYPER/", sampleName,".HYPER.LESIONS.bed", sep = "")))
  {
    lesionAnnotatedSortedHyper <- read.csv2(paste(resultFolder, "/LESIONS_HYPER/", sampleName,".HYPER.LESIONS.bed", sep = ""), sep = "\t", header = FALSE)
  }
  else
  {
    lesionAnnotatedSortedHyper <-data.frame("CHR"="","START"="", "END"="")
    lesionAnnotatedSortedHyper <- lesionAnnotatedSortedHyper[-1,]
  }

  if (file.exists(paste(resultFolder, "/LESIONS_HYPO/", sampleName,".HYPO.LESIONS.bed", sep = "")))
  {
    lesionAnnotatedSortedHypo <- read.csv2(paste(resultFolder, "/LESIONS_HYPO/", sampleName,".HYPO.LESIONS.bed", sep = ""), sep = "\t", header = FALSE)
  }
  else
  {
    lesionAnnotatedSortedHypo <-data.frame("CHR"="","START"="", "END"="")
    lesionAnnotatedSortedHypo <- lesionAnnotatedSortedHypo[-1,]
  }

  lesionAnnotated <- rbind(lesionAnnotatedSortedHyper, lesionAnnotatedSortedHypo)
  colnames(lesionAnnotated) <- c("CHR", "START", "END")
  lesionAnnotatedSorted <- sortByCHRandSTART(lesionAnnotated)
  dumpSampleAsBedFile(
    dataToDump = lesionAnnotatedSorted,
    fileExtension = paste(".", subFileExtension, ".LESIONS.bed", sep = ""),
    resultFolder = resultFolder,
    resultSubFolder = paste("LESIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR", "START", "END", "SAMPLEID")
  )


  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sampleName, " ", "Completed sample ", time_taken)
}
