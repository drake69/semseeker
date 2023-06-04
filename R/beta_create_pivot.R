#' @importFrom doRNG %dorng%
beta_create_pivot <-  function( figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel ) {

  ssEnv <- get_session_info()
  i <- 0
  k <- 0

  keys <- subset(ssEnv$keys,anomalies=="BETA")
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")

  toExport <- c("ssEnv", "tempPopData", "subGroupLabel", "POPULATION", "reportFolder", "mainGroupLabel","sheetList","dir_check_and_create")
  foreach::foreach(k=1:nrow(keys), .export = toExport) %dorng%
    # for(k in 1:nrow(keys))
    {
      grp <- as.character(keys[k,"groups"])
      anomaly <- as.character(keys[k,"anomalies"])
      pivot_file_name <- keys$future_shee_list[k]
      temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
      if(!plyr::empty(temp))
      {
        anomaly <- as.character(keys[k,"anomalies"])
        tempAnomaly <- subset(temp, temp$ANOMALY == as.character(anomaly))
        if(!plyr::empty(tempAnomaly))
        {
          figure <-"BETA"
          anomaly <- "MEAN"
          tempDataFrame <- subset(tempAnomaly, tempAnomaly$FIGURE == figure)
          if(!plyr::empty(tempDataFrame))
          {
            tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula = SAMPLEID + POPULATION ~ KEY, value.var = "VALUE",fun.aggregate = mean, drop = TRUE)

            pivot_subfolder <- dir_check_and_create(reportFolder, anomaly)
            fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
            utils::write.table(t(tempDataFrame), fileName, row.names = T, col.names = F, sep=";")
            tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
            colnames(tempDataFrame) <- tempDataFrame[1,]

            sheet_name <- gsub(" ","", paste0( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
            temp_list <- list(tempDataFrame)

            stats::setNames(temp_list, sheet_name)
          }
        }
      }
    }


}
