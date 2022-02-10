createExcelPivot <-  function(populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel ) {

    HYPO <- NULL
    HYPER <- NULL
    POPULATION <- NULL


    finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)

    # browser()
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
        if(nrow(temp)==0)
          next

        anomaly <- keys[i,2]
        tempAnomaly <- subset(temp, ANOMALY == anomaly )
        if(dim(tempAnomaly)[1]==0)
          next

        figure <- keys[i,3]
        tempDataFrame <- subset(tempAnomaly, FIGURE == figure)
        if(nrow(tempDataFrame)!=0)
        {
          tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + POPULATION ~ KEY, value.var = "freq", sum)
          fileName <- paste0(reportFolder,"/",anomaly,"_",figure, "_", mainGroupLabel,"_", grp,".csv" , sep="")
          write.csv2(t(tempDataFrame), fileName)
          tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
          colnames(tempDataFrame) <- tempDataFrame[1,]

          #store in the first row the names of the sheet
          tempDataFrame[1,] <- gsub(" ","", paste0( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
          list(tempDataFrame)
        }
      }

    # browser()
    for(i in 1:length(sheetList))
    {
      sheetListNames[[i]] <- sheetList[[i]] [1,1]
      sheetList [[i]] <- sheetList[[i]][-1,]
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
