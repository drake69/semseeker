#' sort the pivot table per SAMPLENAME
#'
#' @param pivotTable data to be sorted
#'
#' @return pivt table sorted
#' @export
#'
sortPivot <- function(pivotTable) {

  return(pivotTable)

  # TODO: check the pivot has a column SAMPLENAME
  chr <- colnames(pivotTable)[2:length(colnames(pivotTable))]
  chr <- (gtools::mixedsort(chr))
  pivotTable <- data.frame(pivotTable[, 1], pivotTable[, gtools::mixedsort(chr)])
  colnames(pivotTable)[1] <- "SAMPLENAME"
  pivotTable <- pivotTable[with(pivotTable, gtools::mixedsort(SAMPLENAME)), ]
  return(pivotTable)
}
