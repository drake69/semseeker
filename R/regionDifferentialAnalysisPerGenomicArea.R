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
#' #' @export
#' #'
#' #' @examples
#' #' @importFrom dplyr %>%
#' regionDifferentialAnalysisPerGenomicArea <- function(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel) {
#'   POPULATION <- NULL
#'   i <- NULL
#'   FIGURE <- NULL
#'
#'   chartFolder <- paste(resultFolder, "/Charts/", sep = "")
#'   if (chartFolder != "" && !dir.exists(chartFolder)) {
#'     dir.create(chartFolder)
#'   }
#'
#'   chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel, "/", sep = "")
#'   if (chartFolder != "" && !dir.exists(chartFolder)) {
#'     dir.create(chartFolder)
#'   }
#'
#'   chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel, "/SingleArea", sep = "")
#'   if (chartFolder != "" && !dir.exists(chartFolder)) {
#'     dir.create(chartFolder)
#'   }
#'
#'   logFolder <- paste(resultFolder, "/logs/", sep = "")
#'   if (logFolder != "" && !dir.exists(logFolder)) {
#'     dir.create(logFolder)
#'   }
#'   outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
#'   computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
#'   doParallel::registerDoParallel(computation_cluster)
#'
#'   # options(digits = 22)
#'   parallel::clusterExport(computation_cluster, list("analyzeSingleSample", "dumpSampleAsBedFile", "deltaSingleSample", "createPivotResultFromMultipleBed", "sortByCHRandSTART", "test_match_order", "getLesions", "addCellToDataFrame"))
#'
#'   finalBed <- annotateBed(populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel, resultFolder)
#'
#'   if (is.null(finalBed))
#'     return()
#'
#'   selection_column <- "POPULATION"
#'   finalBed[, subGroupLabel] <- as.factor(finalBed[, subGroupLabel])
#'   finalBed[, mainGroupLabel] <- as.factor(finalBed[, mainGroupLabel])
#'   finalBed$POPULATION <- as.factor(finalBed$POPULATION)
#'   finalBed$FIGURE <- as.factor(finalBed$FIGURE)
#'
#'   subGroups <- levels(as.factor(finalBed[, subGroupLabel]))
#'   mainGroups <- levels(as.factor(finalBed[, mainGroupLabel]))
#'   figures <- levels(as.factor(finalBed$FIGURE))
#'   populations <- c("Reference")
#'
#'   # levels(as.factor(finalBed$POPULATION)) quanti campioni hanno il mainGroup con questo numero di lesioni raggruppati per popolazione
#'
#'   testResult <- data.frame(mainGroupLabel = "", FIGURE = "", subGroupLabel = "", P.VALUE = "", P.VALUE.CORRECTED = "", EXECUTED_TEST = "", POPULATIONS = "")
#'   colnames(testResult) <- c( "FIGURE", subGroupLabel, "P.VALUE", "P.VALUE.CORRECTED", "EXECUTED_TEST", "POPULATIONS")
#'
#'   system(paste0("echo '", paste(colnames(testResult), collapse = "\t"), "' > ", paste(resultFolder, "/Differential_analysis_Area_", mainGroupLabel, ".csv", sep = "")))
#'
#'   pop <- "Reference"
#'   filteredData <- data.frame(subset(finalBed, POPULATION != pop))
#'   keys <- unique(data.frame(subGroupLabel = filteredData[, subGroupLabel], FIGURE = filteredData$FIGURE))
#'   colnames(keys) <- c( subGroupLabel, "FIGURE")
#'   numberOfTest <- dim(keys)[1]
#'
#'   headFile <- paste(c("FIGURE", subGroupLabel, "MEDIAN", "VARIANCE", "SAMPLE_COUNT"), collapse = "\t")
#'   system(paste0("echo '", headFile, "'> ", paste(resultFolder, "/Differential_analysis_Area_", mainGroupLabel, "_Alone.csv", sep = "")))
#'
#'   headFile <- paste(c( "FIGURE", subGroupLabel, "P_VALUE", "P_VALUE_CORRECTED", "EXECUTED_TEST", "POPULATIONS"), collapse = "\t")
#'   system(paste0("echo '", headFile, "' > ", paste(resultFolder, "/Differential_analysis_Area_", mainGroupLabel, ".csv", sep = "")))
#'
#'   headFile <- paste(c( "FIGURE", subGroupLabel, "CASE_MEDIAN", "CASE_VARIANCE", "CONTROL_MEDIAN", "CONTROL_VARIANCE"), collapse = "\t")
#'   system(paste0("echo '", headFile, "' > ", paste(resultFolder, "/Differential_analysis_No_Test_", mainGroupLabel, ".csv", sep = "")))
#'
#'   # foreach::foreach(i = 1:numberOfTest) %dopar% {
#'   for (i in 1:numberOfTest) {
#'     key <- keys[i, ]
#'     subGroup <- key[1, subGroupLabel]
#'     fig <- key$FIGURE
#'
#'     tempData <- data.frame(subset(filteredData, FIGURE == fig))
#'     tempData <- data.frame(subset(tempData, tempData[, subGroupLabel] == subGroup))
#'
#'     fig <- as.character(unique(tempData$FIGURE))
#'     subGroup <- as.character(unique(tempData[1, subGroupLabel]))
#'
#'     freqData <- tempData[, c("freq")]
#'     popData <- tempData[, c("POPULATION")]
#'
#'     if (length(unique(popData)) < 2) {
#'       testResult <- data.frame(matrix(ncol = 6, nrow = 0))
#'       colnames(testResult) <- c( "FIGURE", subGroupLabel, "MEDIAN", "VARIANCE", "SAMPLE_COUNT")
#'       testResult$FIGURE <- fig
#'       testResult[1, subGroupLabel] <- subGroup
#'       testResult[1, "MEDIAN"] <- stats::median(as.vector(tempData$freq))
#'       testResult[1, "VARIANCE"] <- stats::var(as.vector(tempData$freq))
#'       testResult[1, "SAMPLE_COUNT"] <- length(as.vector(tempData$freq))
#'       testResultToDump <- paste(testResult[1, ], collapse = "\t")
#'       system(paste0("echo '", testResultToDump, "' >> ", paste(resultFolder, "/Differential_analysis_Area_", mainGroupLabel, "_Alone.csv", sep = "")))
#'       combination <- paste(fig, "_", paste(unique(tempData$POPULATION), collapse = "_Vs_"), "_", subGroup, sep = "")
#'       print(combination)
#'       next
#'     }
#'
#'
#'     skipCicle <<- FALSE
#'     reswt = tryCatch({
#'       stats::wilcox.test(unlist(freqData) ~ unlist(popData), exact = FALSE)
#'     }, error = function(e) {
#'       print("ERROR")
#'       print(e)
#'       skipCicle <- TRUE
#'     }, finally = {
#'
#'     })
#'
#'     if (skipCicle)
#'       next
#'
#'     if (is.na(reswt$p.value)) {
#'       testResult <- data.frame(matrix(ncol = 7, nrow = 0))
#'       colnames(testResult) <- c(mainGroupLabel, "FIGURE", subGroupLabel, "HYPO_MEDIAN", "HYPO_VARIANCE", "HYPER_MEDIAN", "HYPER_VARIANCE")
#'
#'
#'       testResult$FIGURE <- fig
#'       testResult[1, subGroupLabel] <- subGroup
#'       testResult$CASE_MEDIAN <- stats::median(subset(tempData, POPULATION == "Case")$freq)
#'       testResult$CASE_VARIANCE <- stats::var(subset(tempData, POPULATION == "Case")$freq)
#'       testResult$CONTROL_MEDIAN <- stats::median(subset(tempData, POPULATION == "Control")$freq)
#'       testResult$CONTROL_VARIANCE <- stats::var(subset(tempData, POPULATION == "Control")$freq)
#'       testResultToDump <- paste(testResult[1, ], collapse = "\t")
#'       system(paste0("echo '", testResultToDump, "' >> ", paste(resultFolder, "/Differential_analysis_No_Test_Area_", mainGroupLabel, ".csv", sep = "")))
#'       next
#'     }
#'
#'     testResult <- data.frame(matrix(ncol = 7, nrow = 0))
#'     colnames(testResult) <- c("FIGURE", subGroupLabel, "P_VALUE", "P_VALUE_CORRECTED", "EXECUTED_TEST", "POPULATIONS")
#'     testResult[1,"FIGURE"] <- fig
#'     testResult[1, subGroupLabel] <- subGroup
#'     testResult$EXECUTED_TEST <- numberOfTest
#'     testResult$POPULATIONS <- 2
#'     combination <- paste(fig, "_", paste(unique(tempData$POPULATION), collapse = "_Vs_"), "_", subGroup, sep = "")
#'     print(combination)
#'
#'     testResult$P_VALUE <- reswt$p.value
#'     testResult$P_VALUE_CORRECTED <- 0.05/numberOfTest
#'
#'     if (reswt$p.value < 0.05/numberOfTest) {
#'       testResultToDump <- paste(testResult[1, ], collapse = "\t")
#'       system(paste0("echo '", testResultToDump, "' >> ", paste(resultFolder, "/Differential_analysis_Area_", mainGroupLabel, ".csv", sep = "")))
#'
#'       chart_dataset <- tempData
#'       Groups <- as.factor(chart_dataset[, selection_column])
#'       chart_object <- vector(mode = "list")
#'
#'       chart_object[[i]] <- ggplot2::qplot(
#'         Groups,
#'         chart_dataset[, "freq"],
#'         geom = c("boxplot"),
#'         fill = Groups, alpha = I(0.5),
#'         ylab = paste(fig, " COUNTS")) + ggplot2::ggtitle(paste0(fig, " ", subGroup, sep = "")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_jitter(position = ggplot2::position_jitter(0.2))
#'
#'       ggplot2::ggsave(filename = paste(chartFolder, "/", paste(unique(tempData$POPULATION), collapse = "_Vs_"), "_", fig, "_Box_Plot_", subGroup, ".jpg", sep = ""),
#'                       plot = chart_object[[i]],
#'                       device = NULL,
#'                       path = NULL,
#'                       scale = 1,
#'                       width = NA,
#'                       height = NA,
#'                       units = c("in", "cm", "mm"),
#'                       dpi = 300,
#'                       limitsize = TRUE)
#'     }
#'   }
#'
#'
#' }
