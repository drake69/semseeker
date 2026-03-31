#' #' @importFrom doRNG %dorng%
#' beta_create_pivot <-  function( figures, markers, subGroups, probes_prefix, mainGroupLabel, subGroupLabel ) {
#'
#'   ssEnv <- get_session_info()
#'   i <- 0
#'   k <- 0
#'
#'   keys <- subset(ssEnv$keys,markers=="BETA")
#'   reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
#'
#'   toExport <- c("ssEnv", "tempPopData", "subGroupLabel", "SAMPLE_GROUP", "reportFolder", "mainGroupLabel","sheetList","dir_check_and_create")
#'   foreach::foreach(k=1:nrow(keys), .export = toExport) %dorng%
#'     # for(k in 1:nrow(keys))
#'     {
#'       grp <- as.character(keys[k,"subareas"])
#'       marker <- as.character(keys[k,"markers"])
#'       pivot_file_name <- keys$future_shee_list[k]
#'       temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
#'       if(!plyr::empty(temp))
#'       {
#'         marker <- as.character(keys[k,"markers"])
#'         tempAnomaly <- subset(temp, temp$MARKER == as.character(marker))
#'         if(!plyr::empty(tempAnomaly))
#'         {
#'           figure <-"BETA"
#'           marker <- "MEAN"
#'           tempDataFrame <- subset(tempAnomaly, tempAnomaly$FIGURE == figure)
#'           if(!plyr::empty(tempDataFrame))
#'           {
#'             tempDataFrame <- reshape2::dcast(data = tempDataFrame, formula = SAMPLEID + SAMPLE_GROUP ~ KEY, value.var = "VALUE",fun.aggregate = mean, drop = TRUE)
#'
#'             pivot_subfolder <- dir_check_and_create(reportFolder, marker)
#'             fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
#'             utils::write.table(t(tempDataFrame), fileName, row.names = T, col.names = F, sep=";")
#'             tempDataFrame <- as.data.frame( cbind(colnames(tempDataFrame), t(tempDataFrame)))
#'             colnames(tempDataFrame) <- tempDataFrame[1,]
#'
#'             sheet_name <- gsub(" ","", paste0( marker,"_",figure,"_", mainGroupLabel,"_", grp, sep=""), fixed=TRUE)
#'             temp_list <- list(tempDataFrame)
#'
#'             stats::setNames(temp_list, sheet_name)
#'           }
#'         }
#'       }
#'     }
#'
#'
#' }
