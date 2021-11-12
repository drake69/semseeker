createExcelPivot <-  function(populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel ) {

    HYPO <- NULL
    HYPER <- NULL
    POPULATION <- NULL


    finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)

    browser()
    if (is.null(finalBed))
      return()

    reportFolder <- dir_check_and_create(resultFolderData,"Pivots")

    finalBed <- data.frame(finalBed, "KEY" = finalBed[,mainGroupLabel])
    finalBed[,mainGroupLabel] <- as.factor(finalBed[,mainGroupLabel])
    finalBed[,"KEY"] <- as.factor(finalBed[,"KEY"])
    finalBed[,"FIGURE"] <- as.factor(finalBed[,"FIGURE"])
    finalBed[,"POPULATION"] <- as.factor(finalBed[,"POPULATION"])

    numberOfCase <- length(unique(subset(finalBed, POPULATION == "Case" )$SAMPLEID))
    numberOfControl <- length(unique(subset(finalBed, POPULATION == "Control" )$SAMPLEID))

     # pop <- "Reference"
     # tempPopData <- subset(finalBed, finalBed[,"POPULATION"] != pop)
    tempPopData <- subset(finalBed, freq >0 )
    sheetList <- vector(mode="list")
    sheetListNames <- vector(mode="list")

    # options(digits = 22)
    parallel::clusterExport(envir=environment(), cl = computationCluster, varlist = list(  "sheetList", "sheetListNames"))

    keys <- expand.grid(groups= unique(tempPopData[,subGroupLabel]), "anomalies"= anomalies, "figures"=figures)
    sheetList <- foreach::foreach(i=1:nrow(keys), .export = c("sheetList"), .combine='c', .multicombine=TRUE ) %dopar%
    # for(i in 1:nrow(keys))
      {
        # i <- 1
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

        tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + POPULATION ~ KEY, value.var = "freq", sum)
        fileName <- paste0(reportFolder,"/",anomaly,"_",figure, "_", mainGroupLabel,"_", grp,".csv" , sep="")
        write.csv2(t(tempDataFrame), fileName)
        tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
        colnames(tempDataFrame) <- tempDataFrame[1,]
        list(tempDataFrame[-1,])
      }

    # browser()
    for(i in 1:nrow(keys))
    {
      grp <- keys[i,1]
      anomaly <- keys[i,2]
      figure <- keys[i,3]
      sheetListNames[i] <- gsub(" ","", paste0( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
    }

    fileName <- paste0(reportFolder,"/", mainGroupLabel,".xlsx" , sep="")
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
