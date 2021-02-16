#' takes a bed and its location (build with the details of popuationa nd genomic area)
#' and annoate with detail about genomic area
#' @param populations vector of population to cycle with to build the folder path
#' @param figures vector of hyper /hypo to use to build the folder path
#' @param anomalies vector of lesions/mutations to use to build the folder path
#' @param groups vector of genomic are to cycle and group the annotated data
#' @param probesPrefix prefix to use to get the annotated probes dataset
#' @param columnLabel label of the column of the genomic area gene, island ,dmr etc..
#' @param groupingColumnLabel label of the column of the genomic sub area body, tss1500
#' @param resultFolder folder to which build the folder tree and save the annotated bed
#'
#' @return original bed with genomic area infos
#' @export
#'
#' @examples
#' probesPrefix <- "PROBES_Island_"
#' subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island")
#' mainGroupLabel <- "ISLAND"
#' subGroupLabel <- "RELATION_TO_CPGISLAND"
#' islandBed <- annotateBed(populations ,figures ,anomalies ,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel,resultFolder  )

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

