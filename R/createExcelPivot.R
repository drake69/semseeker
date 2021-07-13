createExcelPivot <-
  function(logFolder, resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel ) {

    HYPO <- NULL
    HYPER <- NULL
    POPULATION <- NULL


    finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel, resultFolder)

    if (is.null(finalBed))
      return()

    reportFolder <- paste(resultFolder, "/Pivots/", sep="")
    if (reportFolder != "" && !dir.exists(reportFolder)) {
      dir.create(reportFolder)
    }

    finalBed <- data.frame(finalBed, "KEY" = finalBed[,mainGroupLabel])
    finalBed[,mainGroupLabel] <- as.factor(finalBed[,mainGroupLabel])
    finalBed[,"KEY"] <- as.factor(finalBed[,"KEY"])
    finalBed[,"FIGURE"] <- as.factor(finalBed[,"FIGURE"])
    finalBed[,"POPULATION"] <- as.factor(finalBed[,"POPULATION"])


    numberOfCase <- length(unique(subset(finalBed, POPULATION == "Case" )$SAMPLENAME))
    numberOfControl <- length(unique(subset(finalBed, POPULATION == "Control" )$SAMPLENAME))

    # pop <- "Reference"
    # tempPopData <- subset(finalBed, finalBed[,"POPULATION"] != pop)
    tempPopData <- subset(finalBed, freq >0 )
    sheetList <- vector(mode="list")
    sheetListNames <- vector(mode="list")
    # i <- 1
    # for (grp in unique(tempPopData[,subGroupLabel]))
    # {
    #   temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
    #   if(dim(temp)[1]==0)
    #     next
    #   for (anomaly in anomalies)
    #   {
    #     tempAnomaly <- subset(temp, ANOMALY == anomaly )
    #     if(dim(tempAnomaly)[1]==0)
    #       next
    #     for (figure in figures)
    #     {
    #       tempDataFrame <- subset(tempAnomaly, FIGURE == figure)
    #       if(dim(tempDataFrame)[1]==0)
    #         next
    #       tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLENAME + POPULATION ~ KEY, value.var = "freq", sum)
    #       # browser()
    #       # row.names(tempDataFrame) <- tempDataFrame$SAMPLENAME
    #       fileName <- paste(reportFolder,"/",anomaly,"_",figure, "_", mainGroupLabel,"_", grp,".csv" , sep="")
    #       write.csv2(t(tempDataFrame), fileName)
    #       sheetList[[i]] <- tempDataFrame
    #       sheetListNames[i] <- gsub(" ","", paste( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
    #       i <- i +1
    #     }
    #   }
    # }
#
#     nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
#     outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
#     print(outFile)
#     computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
#     doParallel::registerDoParallel(computation_cluster)
#
#     # options(digits = 22)
#     parallel::clusterExport(envir=environment(), cl = computation_cluster,
#                             varlist = list(  "sheetList", "sheetListNames"))

    keys <- expand.grid(groups= unique(tempPopData[,subGroupLabel]), "anomalies"= anomalies, "figures"=figures)
    # foreach::foreach(i=1:nrow(keys), .export = c("sheetList", "sheetListNames")) %dopar%
    for(i in 1:nrow(keys))
      {
        grp <- keys[i,1]
        temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
        if(dim(temp)[1]==0)
          next

        anomaly <- keys[i,2]
        tempAnomaly <- subset(temp, ANOMALY == anomaly )
        if(dim(tempAnomaly)[1]==0)
          next

        figure <- keys[i,3]
        tempDataFrame <- subset(tempAnomaly, FIGURE == figure)
        if(dim(tempDataFrame)[1]==0)
          next

        tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLENAME + POPULATION ~ KEY, value.var = "freq", sum)
        fileName <- paste(reportFolder,"/",anomaly,"_",figure, "_", mainGroupLabel,"_", grp,".csv" , sep="")
        write.csv2(t(tempDataFrame), fileName)
        sheetList[[i]] <- tempDataFrame
        sheetListNames[i] <- gsub(" ","", paste( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
      }



    fileName <- paste(reportFolder,"/", mainGroupLabel,".xlsx" , sep="")
    names(sheetList) <- as.vector(sheetListNames)
    try(
      openxlsx::write.xlsx(
        x = sheetList,
        file = fileName,
        asTable = TRUE,
        overwrite = TRUE
      )
    )
    message("Saved spreadsheet file:")
    message(fileName)
}
