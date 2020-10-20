
#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
sortByCHRandSTART <- function(dataframe) {

  dataframe$CHR <- droplevels(dataframe$CHR)

  dataframe$CHR <- factor(dataframe$CHR, levels = mixedsort(attributes(dataframe$CHR)$levels))
  dataframe <- dataframe[with(dataframe, order(CHR, START)), ]
  return(dataframe)
}
