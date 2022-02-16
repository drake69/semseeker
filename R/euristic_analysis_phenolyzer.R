euristic_analysis_phenolizer <- function(resultFolder, terms)
{

   init_env(resultFolder)

  geneAnnotatedFile <- utils::read.csv(file.path(ssEnv$resultFolderData , "GENE_ANNOTATED.csv"))
  geneAnnotatedFile <-subset(geneAnnotatedFile,geneAnnotatedFile$POPULATION != "Reference")
  # semen;azoospermia,sperm;hypertension;thyroid hormones;semen quality

  anomalies <- c("MUTATIONS","LESIONS")
  for (anomaly in anomalies)
  {
    # anomaly <- "LESIONS"
    geneAnnotated <-subset(geneAnnotatedFile,geneAnnotatedFile$ANOMALY == anomaly)

    # tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, SAMPLEID  ~ GENE, value.var = "freq", sum))
    # tempDataFrame <- as.data.frame(t(tempDataFrame))
    # tempDataFrame <- apply(tempDataFrame[],1, as.numeric)
    # cutOff <- mean(tempDataFrame, na.rm = TRUE)
    # geneAnnotated <- geneAnnotated[geneAnnotated$freq>cutOff, ]

    geneAnnotated <-unique(geneAnnotated[, c("GENE", "POPULATION", "ANOMALY", "SAMPLEID")])
    geneAnnotated$freq <- 1

    tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, POPULATION  ~ GENE, value.var = "freq", sum))

    geneMutated <- as.data.frame(t(tempDataFrame))
    names(geneMutated) <- geneMutated[1, ]
    geneMutated <- geneMutated[-1, ]
    geneMutated$gene <- rownames(geneMutated)
    # head(geneMutated[order(geneMutated$Case, decreasing = TRUE),])

    # caseCutOff <- mean(as.numeric(geneMutated$Case))
    # geneMutated <- geneMutated[(as.numeric(geneMutated$Case) > caseCutOff & as.numeric(geneMutated$Control) == 0), ]
    geneMutated <- geneMutated[(as.numeric(geneMutated$Case) > 0 & as.numeric(geneMutated$Control) == 0), ]

    geneMutated$gene <- rownames(geneMutated)
    # utils::write.table(
    #   row.names(geneMutated),
    #   file.path(ssEnv$resultFolderData, paste0("/Euristic/Phenolyzer/", anomaly,"_genes_only_case.csv",sep="")),
    #   row.names = FALSE,
    #   col.names = FALSE ,
    #   quote = FALSE,
    #   sep = "\t"
    # )
    # utils::write.table(
    #   geneMutated[order(geneMutated$Case, decreasing = TRUE),] ,
    #   file.path(ssEnv$resultFolderData, paste0("/Euristic/Phenolyzer/", anomaly, "_genes_only_case_all_details.csv", sep="")),
    #   row.names = FALSE,
    #   col.names = TRUE ,
    #   quote = FALSE,
    #   sep = ";"
    # )
    #
    for (term in terms)
    {

      prio_genes <- phenolyzer_call(term)
      if(is.null(prio_genes))
        next

      geneMutated <- geneMutated[order(geneMutated$Case, decreasing = TRUE),]
      res_prio <- merge(geneMutated, prio_genes, by.x = "gene", by.y = "Gene")

      if(nrow(res_prio)>0)
      {
        res_prio <- res_prio[order(res_prio$Score, decreasing = TRUE), ]
        fileName <- file_path_build (dir_check_and_create(ssEnv$resultFolderData,c("Euristic","Phenolyzer")), c("prioritized",anomaly,term,"mutated_genes"),"csv")
        utils::write.table(res_prio,fileName,row.names = FALSE,col.names = TRUE ,quote = FALSE,sep =";")
      }
    }
  }
}
