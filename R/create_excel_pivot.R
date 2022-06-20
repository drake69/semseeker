#' @importFrom doRNG %dorng%
create_excel_pivot <-  function(envir, populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel ) {

  final_bed <-  annotate_bed( envir=envir,  populations,figures ,anomalies,subGroups ,probes_prefix ,mainGroupLabel,subGroupLabel)

  if (is.null(final_bed))
    return()

  reportFolder <- dir_check_and_create(envir$result_folderData,"Pivots")

  final_bed <- data.frame(final_bed, "KEY" = final_bed[,mainGroupLabel])
  final_bed[,mainGroupLabel] <- as.factor(final_bed[,mainGroupLabel])
  final_bed[,"KEY"] <- as.factor(final_bed[,"KEY"])
  final_bed[,"FIGURE"] <- as.factor(final_bed[,"FIGURE"])
  final_bed[,"POPULATION"] <- as.factor(final_bed[,"POPULATION"])

  numberOfCase <- length(unique(subset(final_bed, final_bed$POPULATION == "Case" )$SAMPLEID))
  numberOfControl <- length(unique(subset(final_bed, final_bed$POPULATION == "Control" )$SAMPLEID))

  tempPopData <- subset(final_bed, final_bed$VALUE != 0 )
  sheetList <- vector(mode="list")
  sheetListNames <- vector(mode="list")

  envir$keys <- expand.grid("groups"= unique(tempPopData[,subGroupLabel]), "anomalies"= anomalies, "figures"=figures)

  toExport <- c("envir", "tempPopData", "subGroupLabel", "POPULATION", "reportFolder", "mainGroupLabel","sheetList","dir_check_and_create")
  sheetList <- foreach::foreach(k=1:nrow(envir$keys), .export = toExport, .combine= c , .multicombine=TRUE ) %dorng%
    # for(k in 1:nrow(envir$keys))
    {
      # k <- 1
      grp <- envir$keys[k,"groups"]
      anomaly <- envir$keys[k,"anomalies"]
      temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
      if(!plyr::empty(temp))
      {
        anomaly <- envir$keys[k,"anomalies"]
        tempAnomaly <- subset(temp, temp$ANOMALY == as.character(anomaly))
        if(!plyr::empty(tempAnomaly))
        {
          figure <- as.character(envir$keys[k,"figures"])
          tempDataFrame <- subset(tempAnomaly, tempAnomaly$FIGURE == figure)
          if(!plyr::empty(tempDataFrame))
          {
            if(anomaly=="DELTAS")
              tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula = SAMPLEID + POPULATION ~ KEY, value.var = "VALUE",
                                               fun.aggregate = mean, drop = TRUE)
            else
              tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula =  SAMPLEID + POPULATION ~ KEY, value.var = "VALUE",
                                               fun.aggregate = sum, drop = TRUE)

            pivot_subfolder <- dir_check_and_create(reportFolder, anomaly)
            fileName <- paste0(pivot_subfolder,"/",anomaly,"_",figure, "_", mainGroupLabel,"_", grp,".csv" , sep="")
            utils::write.csv2(t(tempDataFrame), fileName)
            tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
            colnames(tempDataFrame) <- tempDataFrame[1,]

            sheet_name <- gsub(" ","", paste0( anomaly,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
            temp_list <- list(tempDataFrame)
            setNames(temp_list, sheet_name)
          }
        }
      }
    }

  if(length(sheetList)!=0)
  {

    fileName <- paste0(reportFolder,"/", mainGroupLabel,".xlsx" , sep="")
    try(
      {
        openxlsx::write.xlsx(
          x = sheetList,
          file = fileName,
          asTable = TRUE,
          overwrite = TRUE
        )
        message("Saved spreadsheet file:")
        message(fileName)
      }
    )
  }
  else
  {
    message("No pivot tables to save.")
  }

}
