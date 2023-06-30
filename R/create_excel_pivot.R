#' @importFrom doRNG %dorng%
create_excel_pivot <-  function( sample_groups, figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel ) {

  ssEnv <- get_session_info()
  final_bed <-  annotate_bed(sample_groups =   sample_groups,figures =  figures ,markers =  markers,
                              subareas =   subGroups ,probes_prefix =   probes_prefix , area =  mainGroupLabel, groupingColumnLabel = subGroupLabel)
  i <- 0
  k <- 0
  if (is.null(final_bed) | plyr::empty(final_bed) | nrow(final_bed)==0)
    return()

  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")

  final_bed <- data.frame(final_bed, "KEY" = final_bed[,mainGroupLabel])
  final_bed[,mainGroupLabel] <- as.factor(final_bed[,mainGroupLabel])
  final_bed[,"KEY"] <- as.factor(final_bed[,"KEY"])
  final_bed[,"FIGURE"] <- as.factor(final_bed[,"FIGURE"])
  final_bed[,"SAMPLE_GROUP"] <- as.factor(final_bed[,"SAMPLE_GROUP"])
  final_bed[,"VALUE"] <- as.numeric(final_bed[,"VALUE"])

  numberOfCase <- length(unique(subset(final_bed, final_bed$SAMPLE_GROUP == "Case" )$SAMPLEID))
  numberOfControl <- length(unique(subset(final_bed, final_bed$SAMPLE_GROUP == "Control" )$SAMPLEID))

  tempPopData <- subset(final_bed, final_bed$VALUE != 0 )
  sheetList <- vector(mode="list")
  # sheetListNames <- vector(mode="list")

  fileNameXLS <- paste0(reportFolder,"/", mainGroupLabel,".xlsx" , sep="")
  ssEnv$keys <- unique(tempPopData[,c(subGroupLabel,"MARKER","FIGURE")])
  colnames(ssEnv$keys) <- c("subareas","markers","figures")
  # ssEnv$keys <- expand.grid("subareas"= unique(tempPopData[,subGroupLabel]), "markers"= markers, "figures"=figures)
  ssEnv$keys$future_shee_list <-  unique(paste(final_bed$MARKER,"_",final_bed$FIGURE,"_",  mainGroupLabel,"_",final_bed$AREA, sep=""))

  if(file.exists(fileNameXLS))
  {
    old_sheet_list <- readxl::excel_sheets(fileNameXLS)
    # old_workbook <- openxlsx::loadWorkbook(fileNameXLS)
    # old_sheet_list <- old_workbook$sheet_names
    ssEnv$keys <- ssEnv$keys[!(ssEnv$keys$future_shee_list %in% old_sheet_list), ]
  }


  toExport <- c("ssEnv", "tempPopData", "subGroupLabel", "SAMPLE_GROUP", "reportFolder", "mainGroupLabel","sheetList","dir_check_and_create")
  sheetList <- foreach::foreach(k=1:nrow(ssEnv$keys), .export = toExport, .combine= "c" , .multicombine=TRUE ) %dorng%
  # for(k in 1:nrow(ssEnv$keys))
    {
      # browser()
      # k <- 1
      grp <- as.character(ssEnv$keys[k,"subareas"])
      marker <- as.character(ssEnv$keys[k,"markers"])
      pivot_file_name <- ssEnv$keys$future_shee_list[k]
      temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
      if(!plyr::empty(temp))
      {
        marker <- as.character(ssEnv$keys[k,"markers"])
        tempAnomaly <- subset(temp, temp$MARKER == as.character(marker))
        if(!plyr::empty(tempAnomaly))
        {
          figure <- as.character(ssEnv$keys[k,"figures"])
          if(marker=="BETA")
            figure <- "MEAN"
          tempDataFrame <- subset(tempAnomaly, tempAnomaly$FIGURE == figure)
          if(!plyr::empty(tempDataFrame))
          {
            # if(marker=="DELTAS")
            #   browser()
            if(marker=="BETA")
            {
              tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula = SAMPLEID + SAMPLE_GROUP ~ KEY, value.var = "VALUE",
                fun.aggregate = mean, drop = TRUE)
            }
            else
            tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula =  SAMPLEID + SAMPLE_GROUP ~ KEY, value.var = "VALUE",
                                             fun.aggregate = sum, drop = TRUE)

            pivot_subfolder <- dir_check_and_create(reportFolder, marker)
            fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
            utils::write.table(t(tempDataFrame), fileName, row.names = T, col.names = F, sep=";")
            tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
            colnames(tempDataFrame) <- tempDataFrame[1,]

            sheet_name <- gsub(" ","", paste0( marker,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
            temp_list <- list(tempDataFrame)

            stats::setNames(temp_list, sheet_name)
          }
        }
      }
    }

  if(exists("old_sheet_list") & length(sheetList)!=0)
  {
    old_workbook <- openxlsx::loadWorkbook(fileNameXLS)
    for(t in 1:length(sheetList))
    {
      openxlsx::addWorksheet(old_workbook,names(sheetList[t]))
      openxlsx::writeData(old_workbook,names(sheetList[t]),sheetList[t])
    }
    openxlsx::saveWorkbook(old_workbook,fileNameXLS,overwrite = TRUE)
  } else if(length(sheetList)!=0)
  {
    try(
      {
        openxlsx::write.xlsx(
          x = sheetList,
          file = fileNameXLS,
          asTable = TRUE,
          overwrite = TRUE
        )
        message("INFO: ", Sys.time(), " Saved spreadsheet file:", fileNameXLS)
      }
    )
  }
  else
  {
    if(!exists("old_sheet_list"))
      message("INFO: ", Sys.time(), " No pivot tables to save.")
  }

}
