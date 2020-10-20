#' Calculate stochastic epi mutations from a methylation dataset as outcome report of pivot
#'
#' @param a a value.
#' @param b b value
#' @return files into the result folder with pivot table and bedgraph.
#' @examples
#' fbind(
#' a = 11,
#' b = resultFolderValue
#' )
fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}
