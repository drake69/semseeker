#' createSummaryExcelFromCumulativeBedFile build an excel summary of figures and
#' anomaly staring from the cumulative bed files into the result folder
#'
#' anomaly and figure built tree
#' @param probeFeatures probes information CHR, START, END and rowname
#' as probenames
#' @param sampleSheet sheet of samples with at least Sample_ID column
#'
#' @return a saved excel file

#'
createSummaryExcelFromCumulativeBedFile <-
  function(probeFeatures,
           sampleSheet) {

    lesionsComparison <- NULL
    lesionsComparisonByGene <- NULL

    fileName <- file_path_build(ssEnv$resultFolderData ,"results","xlsx")

    hyperLesionsResult <-
      createPivotResultFromMultipleBed(
        anomalyLabel = "LESIONS",
        figureLable = "HYPER",
        probeFeatures = probeFeatures
      )
    hyperLesions <- hyperLesionsResult$byCHR
    # hyperLesions[hyperLesions == 0] <- ""
    hyperLesions <- sortPivot(pivotTable = hyperLesions)


    hypoLesionsResult <-
      createPivotResultFromMultipleBed(
        anomalyLabel = "LESIONS",
        figureLable = "HYPO",
        probeFeatures = probeFeatures
      )
    hypoLesions <- hypoLesionsResult$byCHR
    # hypoLesions[hypoLesions == 0] <- ""
    hypoLesions <- sortPivot(hypoLesions)

    hyperMutationsResult <-
      createPivotResultFromMultipleBed(
        anomalyLabel = "MUTATIONS",
        figureLable = "HYPER",
        probeFeatures = probeFeatures
      )
    hyperMutations <- hyperMutationsResult$byCHR
    # hyperMutations[hyperMutations == 0] <- ""
    hyperMutations <- sortPivot(hyperMutations)

    hypoMutationsResult <-
      createPivotResultFromMultipleBed(
        anomalyLabel = "MUTATIONS",
        figureLable = "HYPO",
        probeFeatures = probeFeatures
      )
    hypoMutations <- hypoMutationsResult$byCHR
    # hypoMutations[hypoMutations == 0] <- ""
    hypoMutations <- sortPivot(hypoMutations)


    # hypoLesions[is.numeric(hypoLesions)] <- hypoLesions[is.numeric(hypoLesions)] * -1 lesionsComparison <- dplyr::bind_rows(hyperLesions, hypoLesions) lesionsComparison[is.na(lesionsComparison)] <- ''
    # lesionsComparison[lesionsComparison == 0] <- '' lesionsComparison <- sortPivot(lesionsComparison) hyperLesionsResult$byGene[is.numeric(hyperLesionsResult$byGene)] <-
    # hyperLesionsResult$byGene[is.numeric(hyperLesionsResult$byGene)] * -1 lesionsComparisonByGene <- dplyr::bind_rows(hyperLesionsResult$byGene, hypoLesionsResult$byGene) lesionsComparisonByGene[is.na(lesionsComparisonByGene)] <- ''
    # lesionsComparisonByGene[lesionsComparisonByGene == 0] <- ''

    # hyperMutationsResultByGene <- hyperMutationsResult$byGene
    # hyperMutationsResultByGene[hyperMutationsResultByGene == 0] <-
    #   ""
    # hyperLesionsResult$byGene[hyperLesionsResult$byGene == 0] <- ""
    # hypoMutationsResult$byGene[hypoMutationsResult$byGene == 0] <-
    #   ""
    # hypoLesionsResult$byGene[hypoLesionsResult$byGene == 0] <- ""

    sheets <- list(
      SUMMARY = sampleSheet,
      # LESION_COMPARISON = lesionsComparison,
      HYPER_MUTATIONS = hyperMutations,
      HYPER_LESIONS = hyperLesions,
      HYPO_MUTATIONS = hypoMutations,
      HYPO_LESIONS = hypoLesions
      # LESION_COMPARISON_ByGene = lesionsComparisonByGene,
      # HYPER_MUTATIONS_ByGene = hyperMutationsResultByGene,
      # HYPER_LESIONS_ByGene = hyperLesionsResult$byGene,
      # HYPO_MUTATIONS_ByGene = hypoMutationsResult$byGene,
      # HYPO_LESIONS_ByGene = hypoLesionsResult$byGene
    )
    openxlsx::write.xlsx(
      x = sheets,
      file = fileName,
      asTable = TRUE
    )
  }
