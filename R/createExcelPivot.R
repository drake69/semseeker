createExcelPivot <-
  function(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel ) {

    HYPO <- NULL
    HYPER <- NULL
    POPULATION <- NULL

    chartFolder <- paste(resultFolder, "/Charts/", sep="")
    if (chartFolder != "" && !dir.exists(chartFolder)) {
      dir.create(chartFolder)
    }

    chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel,"/", sep="")
    if (chartFolder != "" && !dir.exists(chartFolder)) {
      dir.create(chartFolder)
    }

    chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel,"/Grouped", sep="")
    if (chartFolder != "" && !dir.exists(chartFolder)) {
      dir.create(chartFolder)
    }

    finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel, resultFolder)

    if (is.null(finalBed))
      return()

    finalBed <- data.frame(finalBed, "KEY" = paste(finalBed$POPULATION,"_", finalBed[,mainGroupLabel], sep=""))
    finalBed[,mainGroupLabel] <- as.factor(finalBed[,mainGroupLabel])
    finalBed[,"KEY"] <- as.factor(finalBed[,"KEY"])
    finalBed[,"FIGURE"] <- as.factor(finalBed[,"FIGURE"])
    finalBed[,"POPULATION"] <- as.factor(finalBed[,"POPULATION"])

    numberOfCase <- length(unique(subset(finalBed, POPULATION == "Case" )$SAMPLENAME))
    numberOfControl <- length(unique(subset(finalBed, POPULATION == "Control" )$SAMPLENAME))

    pop <- "Reference"
    tempPopData <- subset(finalBed, finalBed[,"POPULATION"] != pop)
    sheetList <- vectore(mode="list")
    i <- 1
    for (grp in unique(tempPopData[,subGroupLabel]))
    {
      temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
      temp <- reshape2::dcast(data = temp, KEY + POPULATION  ~ FIGURE, value.var = "freq", sum)
      temp <- reshape2::dcast(data = temp, HYPO + HYPER  ~ POPULATION, value.var="HYPER", length)


      for (anomaly in anomalies)
      {

        tempDataFrame <- subset(inputBedDataFrame, ANOMALY == anomaly)
        if(dim(tempDataFrame)[1]==0)
          next
        tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLENAME + POPULATION ~ KEY, value.var = "FREQ", sum)
        row.names(tempDataFrame) <- tempDataFrame$SAMPLENAME
        sheetList[[i]] <- tempDataFrame
        i <- i +1
      }
    }
  openxlsx::write.xlsx(
    x = sheetList,
    file = fileName,
    asTable = TRUE
  )
}
