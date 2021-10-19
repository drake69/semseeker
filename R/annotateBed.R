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

#'
#' @examples
#' probesPrefix <- "PROBES_Island_"
#' subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
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
  resultFolder,
  logFolder)
  {

  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
  outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
  # print(outFile)
  computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computation_cluster)

  # options(digits = 22)
  parallel::clusterExport(envir=environment(), cl = computation_cluster, varlist =c("readMultipleBed","PROBES_Gene_3UTR", "PROBES_Gene_5UTR","PROBES_DMR_DMR","PROBES_Gene_Body",
                                                                                    "PROBES_Gene_TSS1500","PROBES_Gene_TSS200","PROBES_Gene_Whole","PROBES_Gene_ExonBnd","PROBES_Gene_1stExon",
                                                                                    "PROBES_DMR_DMR","PROBES_Island_Island","PROBES_Island_N_Shelf","PROBES_Island_N_Shore","PROBES_Island_Whole"))

  finalBed <- NULL
  bedFileName <- file.path( resultFolder, paste0("/", columnLabel, "_annotatedBed.csv" ,sep=""))

  if(file.exists(bedFileName))
  {
    finalBed <-    utils::read.csv(bedFileName, stringsAsFactors = TRUE)
    return(finalBed)
  }

  keys <-
    expand.grid(
      "POPULATION" = populations,
      "FIGURE" = figures,
      "ANOMALY" = anomalies,
      "GROUP" = groups
    )

  finalBed <- foreach::foreach(g = 1:nrow(keys), .combine = rbind) %dopar%
  {
    anomal <- keys[g,"ANOMALY"]
    pop <- keys[g,"POPULATION"]
    fig <- keys[g,"FIGURE"]
    grp <- keys[g,"GROUP"]

    probes <- get(paste(probesPrefix, grp,sep=""))
    resFolder <- paste0(resultFolder,"/",pop,sep="")
    # tempBed <-
    readMultipleBed( resultFolder = resFolder  , anomalyLabel =  anomal, figureLable =  fig, probeFeatures =  probes, columnLabel =  columnLabel, populationName = pop, groupingColumnLabel= groupingColumnLabel)
    # finalBed <- rbind(finalBed, tempBed)
  }

  # message("Annotated bed:")
  # message(bedFileName)
  utils::write.csv(finalBed,bedFileName, row.names = FALSE)
  return(finalBed)


}

