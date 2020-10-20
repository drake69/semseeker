

#' Title
#'
#' @param pivotTable
#'
#' @return
#' @export
#'
#' @examples
sortPivot <- function(pivotTable) {

  # browser()
  chr <- colnames(pivotTable)[2:length(colnames(pivotTable))]
  chr <- (mixedsort(chr))
  pivotTable <- data.frame(pivotTable[, 1], pivotTable[, mixedsort(chr)])
  colnames(pivotTable)[1] <- "SAMPLENAME"
  pivotTable <- pivotTable[with(pivotTable, mixedsort(SAMPLENAME)), ]
  return(pivotTable)

}
