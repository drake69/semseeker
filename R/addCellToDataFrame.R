#' Find a cell using the row selected with value cellValueSelection in the
#' column colSelection and the column colname
#'
#' @param dataFrame data frame which will be replaced the cell value
#' @param colSelection column containing cell to use for selection
#' @param cellValueSelection value to select the row into colSelection
#' @param colname column containing the cell to be replaced
#' @param cellValue value to replace
#'
#' @return dataframe with the cell replaced
#' @export
#'
addCellToDataFrame <- function(dataFrame,
                               colSelection,
                               cellValueSelection,
                               colname,
                               cellValue) {
  columns_names <- colnames(dataFrame)

  if (!colname %in% columns_names) {
    dataFrame[, colname] <- ""
  }

  row_selector <- dataFrame[, colSelection] == cellValueSelection

  dataFrame[row_selector, colname] <- cellValue
  result <- dataFrame
  return(result)
}
