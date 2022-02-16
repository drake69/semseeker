semseeker_pheno <- function(sampleSheet,
                            methylationData,
                            resultFolder,
                            bonferroniThreshold = 0.05,
                            inferenceDetails = NULL,
                            pheno_term,
                            onlySeed = TRUE) {

  probesFilter <- probes_go_association_phenolizer(pheno_term, onlySeed = TRUE, resultFolder )

  # browser()
  methylationData <- methylationData[ rownames(methylationData) %in% probesFilter, ]

  semseeker(sampleSheet = sampleSheet,
            methylationData = methylationData,
            resultFolder = resultFolder,
            inferenceDetails = inferenceDetails)



}
