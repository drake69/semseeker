#' mutations_get
#'
#' @param values values of methylation
#' @param figure figure to get Mutaions of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison
#' @param probe_features probes features probe, chr, start,end
#' @param sampleName name of the sample
#'
#' @return mutations
#'
#'
mutations_get <- function(values, figure,thresholds, probe_features, sampleName)
{
  ### get probesOverThreshold ################################################################################################

  comparison <- switch(
    figure,
    "HYPO"=`<`,
    "HYPER"=`>`
  )
  # if (!test_match_order(row.names(values), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and Values!", Sys.time())
  # }
  # if (!test_match_order(row.names(thresholds), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and Thresholds!", Sys.time())
  # }

  mutation <- as.numeric(comparison(values, thresholds))


  # names(mutation) <- rownames(values)

  # message(sampleName, " ", "Got probesOverThreshold ", Sys.time())

  ### get mutation_annotated_sorted ################################################################################################

  mutationAnnotated <- data.frame(as.data.frame(probe_features), "MUTATIONS" = mutation, row.names = probe_features$PROBE)

  # if (!test_match_order(row.names(mutationAnnotated), probe_features$PROBE)) {
  #   stop("Wrong order matching Probes and Mutation!", Sys.time())
  # }

  mutation_annotated_sorted <- sort_by_chr_and_start(mutationAnnotated)

  if (!test_match_order(mutation_annotated_sorted$PROBE, row.names(mutation_annotated_sorted))) {
    stop("Mutation annotation sorted is not coherent with probe informations order!")
  }

  result <- c("mutationCount" = sum(mutation_annotated_sorted$MUTATIONS), "lesionCount" = 0, "probesCount" = 0)

  # message(sampleName, " ", "Got mutation_annotated_sorted ", Sys.time())
  return(mutation_annotated_sorted)
}
