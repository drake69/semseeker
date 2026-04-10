#' mutations_get
#'
#' @param values values of methylation
#' @param figure figure to get Mutaions of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison
#' @param sampleName name of the sample
#'
#' @return mutations
#'
#'
mutations_get <- function(values, figure,thresholds, sampleName)
{
  # sort values and thresholds by chr and start
  values <- sort_by_chr_and_start(values)
  thresholds <- sort_by_chr_and_start(thresholds)

  if (figure == "HYPO") {
    mutation <- as.numeric(values[,4] < thresholds$signal_inferior_thresholds)
  }
  if (figure == "HYPER") {
    mutation <- as.numeric(values[,4] > thresholds$signal_superior_thresholds)
  }
  mutationAnnotated <- data.frame("CHR" = thresholds$CHR,"START" = thresholds$START,"END"=thresholds$END, "MUTATIONS" = mutation)
  mutation_annotated_sorted <- sort_by_chr_and_start(mutationAnnotated)

  return(mutation_annotated_sorted)
}
