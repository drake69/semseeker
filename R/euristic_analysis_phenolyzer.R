

euristic_analysis_phenolizer <- function()
{
  geneAnnotated <- read.csv(file.path(resultFolder, "/GENE_annotatedBed.csv"))
  geneAnnotated <-subset(geneAnnotated,geneAnnotated$POPULATION != "Reference")

  tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, SAMPLEID  ~ GENE, value.var = "freq", sum))
  tempDataFrame <- as.data.frame(t(tempDataFrame))
  tempDataFrame <- apply(tempDataFrame[],1, as.numeric)
  cutOff <- mean(tempDataFrame, na.rm = TRUE)
  geneAnnotated <- geneAnnotated[geneAnnotated$freq>cutOff, ]

  geneAnnotated <-unique(geneAnnotated[, c("GENE", "POPULATION", "ANOMALY", "SAMPLEID")])
  geneAnnotated$freq <- 1

  tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, POPULATION  ~ GENE, value.var = "freq", sum))

  geneMutated <- as.data.frame(t(tempDataFrame))
  names(geneMutated) <- geneMutated[1, ]
  geneMutated <- geneMutated[-1, ]
  geneMutated$gene <- rownames(geneMutated)
  # head(geneMutated[order(geneMutated$Case, decreasing = TRUE),])

  caseCutOff <- mean(as.numeric(geneMutated$Case))

  geneMutated <- geneMutated[(as.numeric(geneMutated$Case) > caseCutOff & as.numeric(geneMutated$Control) == 0), ]

  geneMutated$gene <- rownames(geneMutated)
  write.table(
    row.names(geneMutated),
    file.path(resultFolder, "/phenolyzer/mutated_genes_only_case.csv"),
    row.names = FALSE,
    col.names = FALSE ,
    quote = FALSE,
    sep = "\t"
  )
  write.table(
    geneMutated[order(geneMutated$Case, decreasing = TRUE),] ,
    file.path(resultFolder, "/phenolyzer/mutated_genes_only_case_all_details.csv"),
    row.names = FALSE,
    col.names = TRUE ,
    quote = FALSE,
    sep = ";"
  )
  # semen;azoospermia,sperm;hypertension;thyroid hormones;semen quality

  for (term in c(
    "semen",
    "azoospermia",
    "sperm",
    "hypertension",
    "thyroid hormones",
    "semen quality",
    "oligospermia",
    "male reproduction system",
    "reproduction system"
  ))
  {
    # term <- "hypertension"
    #mkdir cd /home/lcorsaro/Desktop/diossina/semseeker/phenolyzer/
    resFolder1 <- file.path(resultFolder, "phenolyzer")

    if (!dir.exists(resFolder))
      dir.create(resFolder1)

    #cd /home/lcorsaro/Desktop/diossina/semseeker/phenolyzer/
    #mkdir hypertension

    resFolder <-file.path(resultFolder, "phenolyzer", gsub(" ", "_", term))

    if (!dir.exists(resFolder))
      dir.create(resFolder)

    filePrioritisedGenes <- file.path(resFolder, paste(term, ".final_gene_list", sep = ""))
    if(!file.exists(filePrioritisedGenes))
    {
      #cd hypertension

      #perl /usr/local/lib/phenolyzer/disease_annotation.pl hypertension -prediction -phenotype -logistic -out hypertension -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15

      # pcommand <- paste("cd ", resFolder ," && perl /usr/local/lib/phenolyzer/disease_annotation.pl '" ,term, "' -p -ph -logistic -out '" ,term, "' -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25", sep="")
      pcommand <-
        paste(
          "cd ",
          resFolder ,
          " && perl /usr/local/lib/phenolyzer/disease_annotation.pl '" ,
          term,
          "' -prediction -phenotype -logistic -out '" ,
          term ,
          "' -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15",
          sep = ""
        )

      system(pcommand)
      }

    prio_genes <-  read.csv(filePrioritisedGenes, sep = "\t")

    geneMutated <- geneMutated[order(geneMutated$Case, decreasing = TRUE),]
    res_prio <- merge(geneMutated, prio_genes, by.x = "gene", by.y = "Gene")

    if(nrow(res_prio)>0)
    {
      res_prio <- res_prio[order(res_prio$Score, decreasing = TRUE), ]
      fileName <- paste("/phenolyzer/prioritized","_",term,"_mutated_genes.csv",sep = "")
      write.table(res_prio,file.path(resultFolder, fileName),row.names = FALSE,col.names = TRUE ,quote = FALSE,sep =";")
    }

  }
}
