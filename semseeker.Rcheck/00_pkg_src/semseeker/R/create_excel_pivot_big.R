#' #' @importFrom doRNG %dorng%
#' create_excel_pivot_big <-  function(ssEnv ) {
#'
#'   k <- 0
#'   reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
#'   sample_sheet <- utils::read.csv2(file.path(ssEnv$result_folderData,"sample_sheet_result.csv"))
#'   sample_names <- sample_sheet$Sample_ID
#'   toExport <- c("ssEnv","dir_check_and_create","sample_names","k")
#'   pivot <- foreach::foreach(k=1:nrow(ssEnv$keys_markers_figures_areas), .export = toExport, .combine= rbind , .multicombine=TRUE ) %dorng%
#'   # for(k in 1:nrow(ssEnv$keys_markers_figures_areas))
#'   {
#'     # k <- 1
#'     marker <- as.character(ssEnv$keys_markers_figures_areas[k,"MARKER"])
#'     figure <- as.character(ssEnv$keys_markers_figures_areas[k,"FIGURE"])
#'     group <- as.character(ssEnv$keys_markers_figures_areas[k,"AREA"])
#'     subgroup <- as.character(ssEnv$keys_markers_figures_areas[k,"SUBAREA"])
#'
#'     pivot_file_name <- paste(marker,"_",figure,"_",  group,"_",subgroup, sep="")
#'     pivot_subfolder <- dir_check_and_create(reportFolder, marker)
#'     fileName <- paste0(pivot_subfolder,"/",pivot_file_name,".csv" , sep="")
#'     if(file.exists(fileName))
#'     {
#'       tempDataFrame <- utils::read.csv2(fileName, header = T)
#'       colnames(tempDataFrame) <- tempDataFrame[1,]
#'       missed_samples <- c("SAMPLEID",sample_names) [ !(c("SAMPLEID",sample_names) %in% colnames(tempDataFrame))]
#'       if(length(missed_samples)>0)
#'       {
#'         # browser()
#'         to_fill <- as.data.frame(matrix(nrow = nrow(tempDataFrame), ncol = length(missed_samples)))
#'         to_fill [is.na(to_fill)] <- 0
#'         colnames(to_fill) <- missed_samples
#'         tempDataFrame <- cbind(tempDataFrame, to_fill)
#'       }
#'       tempDataFrame <- tempDataFrame[-1,c("SAMPLEID",sample_names)]
#'       tempDataFrame <- tempDataFrame[-1,c("SAMPLEID",sample_names)]
#'       colnames(tempDataFrame) <- c("AREA",sample_names)
#'       tempDataFrame$MARKER <- marker
#'       tempDataFrame$FIGURE <- figure
#'       tempDataFrame$AREA <- group
#'       # tempDataFrame$SUBAREA <- subgroup
#'       # if(exists("pivot"))
#'       #   pivot <- rbind(pivot, tempDataFrame)
#'       # else
#'       #   pivot <- tempDataFrame
#'       tempDataFrame
#'     }
#'   }
#'
#'   pivot [is.na(pivot)] <- 0
#'   fileName <- paste0(reportFolder,"/PIVOT.fst" , sep="")
#'   fst::write_fst( x =  pivot ,path =  fileName  , compress = T )
#' }
