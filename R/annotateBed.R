#' Title
#'
#' @param populations
#' @param figures
#' @param anomalies
#' @param groups
#' @param probesPrefix
#' @param columnLabel
#' @param groupingColumnLabel
#' @param resultFolder
#'
#' @return
#' @export
#'
#' @examples
annotateBed <- function (
  populations ,
  figures ,
  anomalies ,
  groups ,
  probesPrefix ,
  columnLabel ,
  groupingColumnLabel,
  resultFolder
  )
  {


  finalBed <- NULL
  bedFileName <- paste0(resultFolder,"/", columnLabel, "_annotatedBed.bed" ,sep="")

  if(file.exists(bedFileName))
  {
    finalBed <-    utils::read.csv(bedFileName, stringsAsFactors = TRUE)
    return(finalBed)
  }

  for (pop in populations) {
    for (fig in figures) {
      for (anomal in anomalies) {
        for (grp in groups)
        {
          probes <- get(paste(probesPrefix, grp,sep=""))
          resFolder <- paste0(resultFolder,"/",pop,sep="")
          tempBed <-  readMultipleBed( resultFolder = resFolder  , anomalyLabel =  anomal, figureLable =  fig, probeFeatures =  probes, columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
          finalBed <- rbind(finalBed, tempBed)


        }
      }
    }
  }

  utils::write.csv(finalBed,bedFileName, row.names = FALSE)
  return(finalBed)


}

