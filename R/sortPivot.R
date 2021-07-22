#' sort the pivot table per SAMPLEID
#'
#' @param pivotTable data to be sorted
#'
#' @return pivt table sorted

#'
sortPivot <- function(pivotTable) {

  return(pivotTable)

  # TODO: check the pivot has a column SAMPLEID
  chr <- colnames(pivotTable)[2:length(colnames(pivotTable))]
  chr <- (gtools::mixedsort(chr))
  pivotTable <- data.frame(pivotTable[, 1], pivotTable[, gtools::mixedsort(chr)])
  colnames(pivotTable)[1] <- "SAMPLEID"
  pivotTable <- pivotTable[with(pivotTable, gtools::mixedsort(SAMPLEID)), ]
  return(pivotTable)
}
