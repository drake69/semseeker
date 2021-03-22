createExcelPivot <- function(inputBedDataFrame, anomalies, groupLabel, groupColumnIDs ,resultFolder) {

  chartFolder <- paste(resultFolder, "/Charts/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  chartFolder <- paste(resultFolder, "/Charts/", groupLabel,"/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  if (is.null(inputBedDataFrame))
    return()

  colnames(inputBedDataFrame) <- c("MAINGROUP","SAMPLENAME", "SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
  inputBedDataFrame$MAINGROUP <- as.factor(inputBedDataFrame$MAINGROUP)
  inputBedDataFrame$SAMPLENAME <- as.factor(inputBedDataFrame$SAMPLENAME)
  inputBedDataFrame$SUBGROUP <- as.factor(inputBedDataFrame$SUBGROUP)
  inputBedDataFrame$FIGURE <- as.factor(inputBedDataFrame$FIGURE)
  inputBedDataFrame$ANOMALY <- as.factor(inputBedDataFrame$ANOMALY)
  inputBedDataFrame$POPULATION <- as.factor(inputBedDataFrame$POPULATION)

  if(length(groupColumnIDs)==2)
  {
    inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame[, groupColumnIDs[2]],"_",inputBedDataFrame$FIGURE,sep=""))
  }
  if(length(groupColumnIDs)==1)
  {
    inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame$FIGURE,sep=""))
  }
  inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)
  inputBedDataFrame <- subset(inputBedDataFrame, POPULATION != "Reference")
  pops <- unique(inputBedDataFrame$POPULATION)

  sheetList <- vectore(mode="list")
  i <- 1
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
  openxlsx::write.xlsx(
    x = sheetList,
    file = fileName,
    asTable = TRUE
  )
}
