#' #' Title
#' #'
#' #' @param resultFolder
#' #' @param populations
#' #' @param figures
#' #' @param anomalies
#' #' @param subGroups
#' #' @param probesPrefix
#' #' @param mainGroupLabel
#' #' @param subGroupLabel
#' #'
#' #' @return
#' 
#' #'
#' #' @examples
#' #' @importFrom dplyr %>%
#' regionDifferentialAnalysis <-
#'   function(resultFolder,
#'            populations,
#'            figures,
#'            anomalies,
#'            subGroups,
#'            probesPrefix,
#'            mainGroupLabel,
#'            subGroupLabel) {
#'     POPULATION <- NULL
#'     i <- NULL
#'     FIGURE <- NULL
#'
#'     chartFolder <- paste0(resultFolder, "/Charts/", sep = "")
#'     if (chartFolder != "" && !dir.exists(chartFolder)) {
#'       dir.create(chartFolder)
#'     }
#'
#'     chartFolder <-
#'       paste0(resultFolder, "/Charts/", mainGroupLabel, "/", sep = "")
#'     if (chartFolder != "" && !dir.exists(chartFolder)) {
#'       dir.create(chartFolder)
#'     }
#'
#'     chartFolder <-
#'       paste0(resultFolder,
#'             "/Charts/",
#'             mainGroupLabel,
#'             "/SingleGene",
#'             sep = "")
#'     if (chartFolder != "" && !dir.exists(chartFolder)) {
#'       dir.create(chartFolder)
#'     }
#'
#'     logFolder <- paste0(resultFolder, "/logs/", sep = "")
#'     if (logFolder != "" && !dir.exists(logFolder)) {
#'       dir.create(logFolder)
#'     }
#'     outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
#'     computationCluster <-
#'       parallel::makeCluster(
#'         parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1,
#'         type = "PSOCK",
#'         outfile = outFile
#'       )
#'     doParallel::registerDoParallel(computationCluster)
#'
#'     # options(digits = 22)
#'     parallel::clusterExport(
#'       computationCluster,
#'       list(
#'         "analyzeSingleSample",
#'         "dumpSampleAsBedFile",
#'         "deltaSingleSample",
#'         "createPivotResultFromMultipleBed",
#'         "sortByCHRandSTART",
#'         "test_match_order",
#'         "getLesions",
#'         "addCellToDataFrame"
#'       )
#'     )
#'
#'     finalBed <-
#'       annotateBed(
#'         populations,
#'         figures ,
#'         anomalies,
#'         subGroups ,
#'         probesPrefix ,
#'         mainGroupLabel,
#'         subGroupLabel,
#'         resultFolder
#'       )
#'
#'     if (is.null(finalBed))
#'       return()
#'
#'     selection_column <- "POPULATION"
#'     finalBed[, subGroupLabel] <- as.factor(finalBed[, subGroupLabel])
#'     finalBed[, mainGroupLabel] <- as.factor(finalBed[, mainGroupLabel])
#'     finalBed$POPULATION <- as.factor(finalBed$POPULATION)
#'     finalBed$FIGURE <- as.factor(finalBed$FIGURE)
#'
#'     subGroups <- levels(as.factor(finalBed[, subGroupLabel]))
#'     mainGroups <-  levels(as.factor(finalBed[, mainGroupLabel]))
#'     figures <- levels(as.factor(finalBed$FIGURE))
#'     populations <- c("Reference")
#'
#'     #levels(as.factor(finalBed$POPULATION))
#'     # quanti campioni hanno il mainGroup con questo numero di lesioni raggruppati per popolazione
#'
#'     testResult <-
#'       data.frame(
#'         mainGroupLabel = "",
#'         "FIGURE" = "",
#'         subGroupLabel = "",
#'         "P.VALUE" = "",
#'         "P.VALUE.CORRECTED" = "",
#'         "EXECUTED_TEST" = "",
#'         "POPULATIONS" = ""
#'       )
#'     colnames(testResult) <-
#'       c(
#'         mainGroupLabel,
#'         "FIGURE",
#'         subGroupLabel,
#'         "P.VALUE",
#'         "P.VALUE.CORRECTED",
#'         "EXECUTED_TEST",
#'         "POPULATIONS"
#'       )
#'
#'     system(paste0(
#'       "echo '",
#'       paste0(colnames(testResult), collapse = "\t"),
#'       "' > ",
#'       paste0(
#'         resultFolder,
#'         "/Differential_analysis_",
#'         mainGroupLabel,
#'         ".csv",
#'         sep = ""
#'       )
#'     ))
#'
#'     pop <- "Reference"
#'     filteredData <- data.frame(subset(finalBed, POPULATION != pop))
#'     keys <-
#'       unique(
#'         data.frame(
#'           mainGroupLabel = filteredData[, mainGroupLabel],
#'           subGroupLabel = filteredData[, subGroupLabel],
#'           "FIGURE" = filteredData$FIGURE
#'         )
#'       )
#'     colnames(keys) <- c(mainGroupLabel, subGroupLabel, "FIGURE")
#'     numberOfTest <- dim(keys)[1]
#'
#'     # fd <- filteredData %>%  dplyr::group_by(filteredData[,mainGroupLabel],filteredData[,subGroupLabel],filteredData$FIGURE)
#'     # fd <- dplyr::group_split(fd)
#'     # filterOnlye <-  parallel::parLapply(computationCluster, fd, function(i) unique(length(i$POPULATION))<2)
#'     # numberOfTest <- length(fd)
#'
#'     # tfun <- function(tempData){
#'     #
#'     #   mainGroup <- unique(tempData[,mainGroupLabel])
#'     #   fig <- unique(tempData$FIGURE)
#'     #   subGroup <- unique(tempData[,subGroupLabel])
#'     #   pop <- "Reference"
#'     #
#'     #   testResult <- data.frame(mainGroupLabel="","FIGURE"="",subGroupLabel="","P.VALUE"="","P.VALUE.CORRECTED"="", "EXECUTED_TEST" ="", "POPULATIONS"="")
#'     #
#'     #   freqData <- tempData[,c("freq")]
#'     #   popData <- tempData[,c("POPULATION")]
#'     #   combination <-   paste0(fig, "_", pop, "_", mainGroup, "_", subGroup, sep="")
#'     #
#'     #   testResult[,mainGroupLabel] <-mainGroup
#'     #   testResult$FIGURE <- fig
#'     #   testResult[,subGroupLabel] <- subGroup
#'     #   testResult$EXECUTED_TEST <- numberOfTest
#'     #   testResult$POPULATIONS <- 2
#'     #
#'     #   print(combination)
#'     #   skipCicle <- FALSE
#'     #   reswt = tryCatch({
#'     #     stats::wilcox.test(unlist(freqData) ~ unlist(popData), exact=FALSE)
#'     #   },error = function(e) {
#'     #     print("ERROR")
#'     #     print(e)
#'     #     skipCicle <- TRUE
#'     #   }, finally = {
#'     #
#'     #   })
#'     #
#'     #   if(skipCicle)
#'     #     next
#'     #
#'     #   testResult$P.VALUE <- reswt$p.value
#'     #   testResult$P.VALUE.CORRECTED <- 0.05/numberOfTest
#'     #
#'     #   if(is.na(reswt$p.value))
#'     #     next
#'     #
#'     #   try(
#'     #     if(reswt$p.value < 0.05/numberOfTest)
#'     #     {
#'     #       system(paste0( "echo '",paste0( testResult[1,], collapse = "\t"), "' >> ", paste0(resultFolder,"/Differential_analysis_Gene.csv", sep = "")))
#'     #
#'     #       chart_dataset <- tempData
#'     #       Groups <- as.factor(chart_dataset[, selection_column])
#'     #       chart_object <- vector(mode = "list")
#'     #
#'     #       chart_object[[i]] <-
#'     #         ggplot2::qplot(
#'     #           Groups,
#'     #           chart_dataset[, "freq"],
#'     #           geom = c("boxplot"),
#'     #           fill = Groups,
#'     #           alpha = I(0.5),
#'     #           ylab = paste0 ( fig, " COUNTS")
#'     #         ) + ggplot2::ggtitle(paste0(fig, " " , mainGroup, sep="")) +  ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_jitter(position =  ggplot2::position_jitter(0.2))
#'     #
#'     #       ggplot2::ggsave(
#'     #         filename = paste0(chartFolder,pop ,"_" ,fig ,  "_Box_Plot_Gene_", mainGroup, "_", subGroup, ".jpg",sep=""),
#'     #         plot = chart_object[[i]],
#'     #         device = NULL,
#'     #         path = NULL,
#'     #         scale = 1,
#'     #         width = NA,
#'     #         height = NA,
#'     #         units = c("in", "cm", "mm"),
#'     #         dpi = 300,
#'     #         limitsize = TRUE
#'     #       )
#'     #     }
#'     #   )
#'     #
#'     # }
#'     # parallel::parLapply(cl = computationCluster, fd, fun = tfun(fd) )
#'
#'     headFile <- paste0(c(mainGroupLabel,"FIGURE",subGroupLabel,"MEDIAN","VARIANCE","SAMPLE_COUNT"),collapse ="\t")
#'     system(paste0("echo '",headFile,"'> ",paste0(resultFolder,"/Differential_analysis_",mainGroupLabel ,"_Alone.csv",sep = "")))
#'
#'     headFile <-paste0( c(mainGroupLabel,"FIGURE",subGroupLabel,"P_VALUE","P_VALUE_CORRECTED","EXECUTED_TEST","POPULATIONS"),collapse="\t")
#'     system(paste0("echo '",headFile,"' > ",paste0(resultFolder,"/Differential_analysis_",mainGroupLabel ,".csv",sep = "")))
#'
#'     headFile <-paste0( c(mainGroupLabel,"FIGURE",subGroupLabel,"CASE_MEDIAN","CASE_VARIANCE","CONTROL_MEDIAN","CONTROL_VARIANCE"),collapse="\t")
#'     system(paste0("echo '",headFile,"' > ",paste0(resultFolder,"/Differential_analysis_No_Test_",mainGroupLabel ,".csv",sep = "")))
#'
#'     # foreach::foreach(i = 1:numberOfTest) %dopar% {
#'     for (i in 1:numberOfTest) {
#'         key <- keys[i, ]
#'         mainGroup <- key[1, mainGroupLabel]
#'         subGroup <- key[1, subGroupLabel]
#'         fig <- key$FIGURE
#'
#'         tempData <- data.frame(subset(filteredData, FIGURE == fig))
#'         tempData <- data.frame(subset(tempData, tempData[, mainGroupLabel] == mainGroup))
#'         tempData <- data.frame(subset(tempData, tempData[, subGroupLabel] == subGroup))
#'
#'         mainGroup <- as.character(tempData[1, mainGroupLabel])
#'         fig <- as.character(unique(tempData$FIGURE))
#'         subGroup <- as.character(unique(tempData[1, subGroupLabel]))
#'
#'         freqData <- tempData[, c("freq")]
#'         popData <- tempData[, c("POPULATION")]
#'
#'         if (length(unique(popData)) < 2)
#'         {
#'           testResult <- data.frame(matrix(ncol = 6, nrow = 0))
#'           colnames(testResult) <-
#'             c(
#'               mainGroupLabel,
#'               "FIGURE",
#'               subGroupLabel,
#'               "MEDIAN",
#'               "VARIANCE",
#'               "SAMPLE_COUNT"
#'             )
#'           testResult[1, mainGroupLabel] <- mainGroup
#'           testResult$FIGURE <- fig
#'           testResult[1, subGroupLabel] <- subGroup
#'           testResult[1,"MEDIAN"] <- stats::median(as.vector(tempData$freq))
#'           testResult[1,"VARIANCE"] <- stats::var(as.vector(tempData$freq))
#'           testResult[1,"SAMPLE_COUNT"] <- length(as.vector(tempData$freq))
#'           testResultToDump <-
#'             paste0(testResult[1, ], collapse = "\t")
#'           system(paste0(
#'             "echo '",
#'             testResultToDump,
#'             "' >> ",
#'             paste0(
#'               resultFolder,
#'               "/Differential_analysis_",
#'               mainGroupLabel ,
#'               "_Alone.csv",
#'               sep = ""
#'             )
#'           ))
#'           combination <-  paste0(fig, "_", paste0(unique(tempData$POPULATION), collapse = "_Vs_"), "_", mainGroup, "_", subGroup, sep = "")
#'           print(combination)
#'           next
#'         }
#'
#'
#'         skipCicle <<- FALSE
#'         reswt = tryCatch({
#'           stats::wilcox.test(unlist(freqData) ~ unlist(popData), exact = FALSE)
#'         }, error = function(e) {
#'           print("ERROR")
#'           print(e)
#'           skipCicle <- TRUE
#'         }, finally = {
#'
#'         })
#'
#'         if (skipCicle)
#'           next
#'
#'         if (is.na(reswt$p.value))
#'           {
#'           testResult <- data.frame(matrix(ncol = 7, nrow = 0))
#'           colnames(testResult) <-
#'             c(
#'               mainGroupLabel,
#'               "FIGURE",
#'               subGroupLabel,
#'               "HYPO_MEDIAN",
#'               "HYPO_VARIANCE",
#'               "HYPER_MEDIAN",
#'               "HYPER_VARIANCE"
#'             )
#'
#'
#'           testResult[1, mainGroupLabel] <- mainGroup
#'           testResult$FIGURE <- fig
#'           testResult[1, subGroupLabel] <- subGroup
#'           testResult$CASE_MEDIAN <- stats::median(subset(tempData, POPULATION =="Case")$freq)
#'           testResult$CASE_VARIANCE <-  stats::var(subset(tempData, POPULATION =="Case")$freq)
#'           testResult$CONTROL_MEDIAN <-  stats::median(subset(tempData, POPULATION =="Control")$freq)
#'           testResult$CONTROL_VARIANCE <-  stats::var(subset(tempData, POPULATION =="Control")$freq)
#'           testResultToDump <-
#'             paste0(testResult[1, ], collapse = "\t")
#'           system(paste0(
#'             "echo '",
#'             testResultToDump,
#'             "' >> ",
#'             paste0(
#'               resultFolder,
#'               "/Differential_analysis_No_Test_",
#'               mainGroupLabel ,
#'               ".csv",
#'               sep = ""
#'             )
#'           ))
#'           next
#'         }
#'
#'         testResult <- data.frame(matrix(ncol = 7, nrow = 0))
#'         colnames(testResult) <-
#'           c(
#'             mainGroupLabel,
#'             "FIGURE",
#'             subGroupLabel,
#'             "P_VALUE",
#'             "P_VALUE_CORRECTED",
#'             "EXECUTED_TEST",
#'             "POPULATIONS"
#'           )
#'         testResult[1, mainGroupLabel] <- mainGroup
#'         testResult$FIGURE <- fig
#'         testResult[1, subGroupLabel] <- subGroup
#'         testResult$EXECUTED_TEST <- numberOfTest
#'         testResult$POPULATIONS <- 2
#'         combination <-  paste0(fig, "_", paste0(unique(tempData$POPULATION), collapse = "_Vs_"), "_", mainGroup, "_", subGroup, sep = "")
#'         print(combination)
#'
#'         testResult$P_VALUE <- reswt$p.value
#'         testResult$P_VALUE_CORRECTED <- 0.05 / numberOfTest
#'
#'         if (reswt$p.value < 0.05 / numberOfTest)
#'         {
#'           testResultToDump <- paste0(testResult[1, ], collapse = "\t")
#'           system(paste0(
#'             "echo '",
#'             testResultToDump,
#'             "' >> ",
#'             paste0(
#'               resultFolder,
#'               "/Differential_analysis_",
#'               mainGroupLabel ,
#'               ".csv",
#'               sep = ""
#'             )
#'           ))
#'
#'           chart_dataset <- tempData
#'           Groups <- as.factor(chart_dataset[, selection_column])
#'           chart_object <- vector(mode = "list")
#'
#'           chart_object[[i]] <-
#'             ggplot2::qplot(
#'               Groups,
#'               chart_dataset[, "freq"],
#'               geom = c("boxplot"),
#'               fill = Groups,
#'               alpha = I(0.5),
#'               ylab = paste0 (fig, " COUNTS")
#'             ) + ggplot2::ggtitle(paste0(fig, " " , mainGroup, sep = "")) +  ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_jitter(position =  ggplot2::position_jitter(0.2))
#'
#'           ggplot2::ggsave(
#'             filename = paste0(
#'               chartFolder,
#'               "/",
#'               paste0(unique(tempData$POPULATION), collapse = "_Vs_"),
#'               "_" ,
#'               fig ,
#'               "_Box_Plot_",
#'               mainGroupLabel ,
#'               "_",
#'               mainGroup,
#'               "_",
#'               subGroup,
#'               ".jpg",
#'               sep = ""
#'             ),
#'             plot = chart_object[[i]],
#'             device = NULL,
#'             path = NULL,
#'             scale = 1,
#'             width = NA,
#'             height = NA,
#'             units = c("in", "cm", "mm"),
#'             dpi = 300,
#'             limitsize = TRUE
#'           )
#'         }
#'     }
#'
#'     #foreach::foreach(i = 1:numberOfTest) %dopar% {
#'     # for (i in 1:numberOfTest) {
#'     #
#'     #     key <- keys[i,]
#'     #     mainGroup <- key[,mainGroupLabel]
#'     #     subGroup <- key[,subGroupLabel]
#'     #     fig <- key$FIGURE
#'     #
#'     #     tempData <- data.frame(subset(filteredData, FIGURE == fig))
#'     #     tempData <- data.frame(subset(tempData, mainGroupLabel == mainGroup))
#'     #     tempData <- data.frame(subset(tempData, subGroupLabel == subGroup))
#'     #     tfun(tempData)
#'     # }
#'     # write.csv(testResults,paste0(resultFolder,"/Differential_analysis_Gene.csv"))
#'
#'   }
