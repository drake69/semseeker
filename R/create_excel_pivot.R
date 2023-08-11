#' @importFrom doRNG %dorng%
create_excel_pivot <-  function() {

  ssEnv <- get_session_info()

  i <- 0
  k <- 0

  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  localKeys <- ssEnv$keys_areas_subareas_markers_figures
  areas <- ssEnv$keys_areas
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))
  else
    progress_bar <- ""

  for (a in 1:length(areas))
  {

    area <- as.character(areas[a])
    tempKeys <- unique(subset(localKeys,localKeys$AREA==area))
    # tempKeys$COMBINED <-  unique(paste(tempKeys$MARKER,"_",tempKeys$FIGURE,"_",tempKeys$AREA,"_",tempKeys$SUBAREA, sep=""))

    sheetList <- vector(mode="list")
    fileNameXLS <- paste0(reportFolder,"/", area,".xlsx" , sep="")

    if(file.exists(fileNameXLS))
    {
      old_sheet_list <- readxl::excel_sheets(fileNameXLS)
      tempKeys <- tempKeys[!(tempKeys$COMBINED %in% old_sheet_list), ]
    }

    toExport <- c("ssEnv", "annotatedData", "subGroupLabel", "SAMPLE_GROUP", "reportFolder", "area",
      "sheetList","dir_check_and_create","tempKeys",
      "progress_bar","progression_index", "progression", "progressor_uuid",
      "owner_session_uuid", "trace","file_path_build", "read_annotated_bed")
    sheetList <- foreach::foreach(k=1:nrow(tempKeys), .export = toExport, .combine= "c" , .multicombine=TRUE ) %dorng%
    # for(k in 1:nrow(tempKeys))
    {
      # browser()
      # k <- 1
      subarea <-  as.character(tempKeys[k,"SUBAREA"])
      marker <- as.character(tempKeys[k,"MARKER"])
      figure <-  as.character(tempKeys[k,"FIGURE"])
      pivot_file_name <- tempKeys$COMBINED[k]

      annotatedData <-  read_annotated_bed(figure,marker,area,subarea)
      annotatedData <- subset(annotatedData, annotatedData$VALUE != 0 )

      # if (is.null(annotatedData) | plyr::empty(annotatedData))
      #  next

      # if (nrow(annotatedData)==0)
      #   next
      #
      if(!plyr::empty(annotatedData))
      {
          if(marker=="BETA")
          {
            annotatedData <- reshape2::dcast(data = annotatedData, formula = SAMPLEID + SAMPLE_GROUP ~ AREA, value.var = "VALUE",
              fun.aggregate = mean, drop = TRUE)
          }
          else
            annotatedData <- reshape2::dcast(data = annotatedData, formula =  SAMPLEID + SAMPLE_GROUP ~ AREA, value.var = "VALUE",
              fun.aggregate = sum, drop = TRUE)

          pivot_subfolder <- dir_check_and_create(reportFolder, marker)
          fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
          utils::write.table(t(annotatedData), fileName, row.names = T, col.names = F, sep=";")
          annotatedData <- as.data.frame( cbind(colnames(annotatedData), t(annotatedData)))
          colnames(annotatedData) <- annotatedData[1,]


          if(ssEnv$showprogress)
            progress_bar(sprintf("Creating pivot table."))

          temp_list <- list(annotatedData)
          names(temp_list) <- c(pivot_file_name)
          if(ssEnv$showprogress)
            progress_bar(sprintf("Creating pivot table."))
          temp_list
          # sheetList <- c(sheetList, temp_list)
      }
      else
        if(ssEnv$showprogress)
          progress_bar(sprintf("Creating pivot table."))

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
    } else
      if(length(sheetList)!=0)
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

}
