#' mutations_get
#'
#' @param values values of methylation
#' @param figure figure to get Mutaions of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison
#' @param probe_features probe_features features probe, chr, start,end
#' @param sampleName name of the sample
#'
#' @return mutations
#'
#'
mutations_get <- function(values, figure,thresholds, probe_features, sampleName)
{
  comparison <- switch(
    figure,
    "HYPO"=`<`,
    "HYPER"=`>`
  )

  mutation <- as.numeric(comparison(values, thresholds))
  mutationAnnotated <- data.frame(as.data.frame(probe_features), "MUTATIONS" = mutation, row.names = probe_features$PROBE)
  mutation_annotated_sorted <- sort_by_chr_and_start(mutationAnnotated)

  if (!test_match_order(mutation_annotated_sorted$PROBE, row.names(mutation_annotated_sorted))) {
    stop("Mutation annotation sorted is not coherent with probe informations order!")
  }
  return(mutation_annotated_sorted)
}
