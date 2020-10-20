#' Find a cell using the row selected with value cellValueSelection in the
#' column colSelection and the column colname
#'
#' @param dataFrame data frame which will be replaced the cell value
#' @param colSelection column containing cell to use for selection
#' @param cellValueSelection value to select the row into colSelection
#' @param colname column containing the cell to be replaced
#' @param cellValue value to replace
#'
#' @return
#' @export
#'
#' @examples
#' addCellToDataFrame (
#' dataFrame = df,
#' colSelection = "Name Surname",
#' cellValueSelection = "Joe Doe",
#' colname = "Age",
#' cellValue = 21
#' )
addCellToDataFrame <- function(dataFrame, colSelection, cellValueSelection, colname, cellValue) {

  if (!colname %in% colnames(dataFrame))
    dataFrame[, colname] <- ""

  dataFrame[dataFrame[, colSelection] == cellValueSelection, colname] <- cellValue

  return(dataFrame)
}

