#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param slidingWindowSize size of the sliding widows to compute epilesions
#' default 11 probes.
#' @param resultFolder folder to save computed epimutations and bedgraphs files.
#' @param logFolder folder to output log
#' @param popStudyDefinition  data frame to select, from the sample sheet,
#' samples to use as control as study population and as refereces two vectors
#' within the first vector the names of the selection colum and tha second
#' vector the study population selector,
#' @param popControlSelector name of samplesheet's column to use as control
#' population selector followed by selection value,
#' @param bonferroniThreshold = 0.05 #threshold to define which pValue
#' accept for lesions definition
#' @return files into the result folder with pivot table and bedgraph.
#' @export
#' @examples
#' calculate(
#' slidingWindowSize = 11,
#' resultFolder = resultFolderValue,
#' logFolder = logFolderValue,
#' popStudyDefinition  = ,
#' popControlSelector = "Case_Ctrl == 0",
#' bonferroniThreshold = 0.05
#' )
calculate <-
    function(sampleSheet,
             methylationData,
             probeFeatures,
             resultFolder,
             logFolder,
             slidingWindowSize = 11,
             popControlSelector,
             popStudyDefinition,
             bonferroniThreshold = 0.05)
        {

        # TODO: check probe feautures has necessary fields CHR, START

        needColumn <- c("Sample_Name")
        missedColumns <-
            needColumn[!(needColumn %in% colnames(sampleSheet))]

        if (length(missedColumns) > 0) {
            stop(
                "File:",
                sampleSheets[1],
                " Lost following columns ",
                missedColumns,
                " ",
                Sys.time()
            )
        }

        if (logFolder != "" && !dir.exists(logFolder)) {
            dir.create(logFolder)
        }

        # reference population
        populations <-
            rbind(
                popStudyDefinition,
                data.frame(LABELS = "CONTROL", SELECTORS = popControlSelector)
            )

        # control population
        controlPopulationSampleSheet <-
            subset(sampleSheet, eval(parse(text = popControlSelector)))
        controlPopulationMatrix <-
            data.frame(PROBE = row.names(methylationData), methylationData[, controlPopulationSampleSheet$Sample_Name])

        if (empty(controlPopulationMatrix) |
            dim(controlPopulationMatrix)[2] < 2) {
            message("Empty methylationData ", Sys.time())
            stop("Empty methylationData ")
        }

        populationControlRangeBetaValues <-
            rangeBetaValuePerProbeAsNormalizedDistribution(controlPopulationMatrix, iqrTimes = 3)


        for (i in 1:dim(populations)[1]) {
            if (is.null(populations[i, "SELECTORS"])) {
                next
            }

            popSelector <- as.character(populations[i, "SELECTORS"])
            popName <- as.character(populations[i, "LABELS"])

            # browser()
            populationSampleSheet <-
                subset(sampleSheet, eval(parse(text = popSelector)))
            populationMatrix <-
                data.frame(PROBE = row.names(methylationData), methylationData[, populationSampleSheet$Sample_Name])

            if (empty(populationMatrix) |
                dim(populationMatrix)[2] < 2) {
                message("Population with Selector ",
                        popSelector,
                        " is empty ",
                        Sys.time())
                next
            }

            analizePopulation(
                populationMatrix = populationMatrix,
                slidingWindowSize = slidingWindowSize,
                resultFolder = paste0(resultFolder, "_", j, sep = ""),
                logFolder = logFolder,
                betaSuperiorThresholds = populationControlRangeBetaValues$betaSuperiorThresholds,
                betaInferiorThresholds = populationControlRangeBetaValues$betaInferiorThresholds,
                sampleSheet = populationSampleSheet,
                probeFeatures = probeFeatures,
                betaMeans = populationControlRangeBetaValues$betaMedianValues,
                populationName = popName,
                bonferroniThreshold = bonferroniThreshold
            )

            rm(populationSampleSheet)
            rm(populationMatrix)
        }

        rm(populationControlRangeBetaValues)
    }
