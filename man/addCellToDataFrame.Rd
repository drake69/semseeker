% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addCellToDataFrame.R
\name{addCellToDataFrame}
\alias{addCellToDataFrame}
\title{Find a cell using the row selected with value cellValueSelection in the
column colSelection and the column colname}
\usage{
addCellToDataFrame(
  dataFrame,
  colSelection,
  cellValueSelection,
  colname,
  cellValue
)
}
\arguments{
\item{dataFrame}{data frame which will be replaced the cell value}

\item{colSelection}{column containing cell to use for selection}

\item{cellValueSelection}{value to select the row into colSelection}

\item{colname}{column containing the cell to be replaced}

\item{cellValue}{value to replace}
}
\value{
dataframe with the cell replaced
}
\description{
Find a cell using the row selected with value cellValueSelection in the
column colSelection and the column colname
}
