#' getMutations
#'
#' @param values values of methylation
#' @param figure figure to get Mutaions of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison
#' @param probeFeatures probes features probe, chr, start,end
#' @param sampleName name of the sample
#'
#' @return mutations
#'
#'
getMutations <- function(values, figure,thresholds, probeFeatures, sampleName)
{
  ### get probesOverThreshold ################################################################################################

  comparison <- switch(
    figure,
    "HYPO"=`<`,
    "HYPER"=`>`
  )

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

  message(sampleName, " ", "Got mutationAnnotatedSorted ", Sys.time())
  return(mutationAnnotatedSorted)
}
