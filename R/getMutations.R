#' getMutations
#'
#' @param values values of methylation
#' @param thresholds threshold to use for comparison
#' @param comparison function to use for compare
#' @param probeFeatures probes features probe, chr, start,end
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
