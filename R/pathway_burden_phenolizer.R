pathway_burden_phenolizer <- function(resultFolder, terms)
{

  init_env(resultFolder)

  figures <- c("BOTH")
  anomalies <- c("MUTATIONS")
  subGroups <- c("Whole")

  probesPrefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"

  geneAnnotatedFile <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)

  geneAnnotatedFile$GENE <- string_normalize(geneAnnotatedFile$GENE)
  geneAnnotatedFile <-subset(geneAnnotatedFile,POPULATION != "Reference")

  anomalies <- c("MUTATIONS")
  # anomalies <- c("MUTATIONS","LESIONS")
  # term <- "spermatogenesis"
  samples <- read.csv(file_path_build(baseFolder =  resultFolderData,detailsFilename = c("sample","sheet","result"),"csv"), sep=";")
  for (anomaly in anomalies)
  {
    # anomaly <- "MUTATIONS"
    geneAnnotated <-subset(geneAnnotatedFile,ANOMALY == anomaly & FIGURE=="BOTH" & GROUP=="Whole")
    for (term in terms)
    {

      prio_genes <- phenolyzer_call(term)
      if(is.null(prio_genes))
        next

      prio_genes$Gene <- string_normalize(prio_genes$Gene)
      # prio_genes <- prio_genes[prio_genes$Status=="SeedGene", ]
      tempDataFrame <- geneAnnotated[ geneAnnotated$GENE %in% prio_genes$Gene, ]
      tempDataFrame$path <- term
      tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID ~ path, value.var = "freq", sum)
      if(!exists("result"))
        result <- merge(samples, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID",all.x=TRUE)
      else
        result <- merge(result, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID",all.x=TRUE)
    }
  }
  fileName <- file_path_build (dir_check_and_create(resultFolder,c("Pathway","Phenolyzer")), c("mutations","both",term,"burden","pathway"),"csv")
  write.table(result,fileName,row.names = FALSE,col.names = TRUE ,quote = FALSE,sep =";")
}
