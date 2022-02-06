#' inference_pathway_analysis_phenolizer calculate and att to sample sheet result
#' the number of mutations or the lesioned gene per each sample
#' @param resultFolder
#' @param terms
#'
#' @return
#'
#' @examples
inference_pathway_analysis_phenolizer <- function(resultFolder,terms)
{

  init_env(resultFolder)

  geneAnnotatedFile <- read.csv(file.path(resultFolderData, "/GENE_ANNOTATED.csv"))
  geneAnnotatedFile <-subset(geneAnnotatedFile,POPULATION != "Reference")

  anomalies <- c("LESIONS","MUTATIONS")
  for(anomaly in anomalies){

    # anomaly <- "LESIONS"
    genesAffected <-subset(geneAnnotatedFile,ANOMALY == anomaly)
    if(anomaly=="LESIONS")
    {
      # one lesions per gene allowed
      genesAffected$freq <-1
      genesAffected <- unique(genesAffected[,c("SAMPLEID","GENE","freq")] )
    }
    genesAffected <- as.data.frame(reshape2::dcast(data = genesAffected, SAMPLEID  ~ GENE, value.var = "freq", sum))

    for (term in terms)
    {
      prio_genes <- phenolyzer_call(term)
      if(is.null(prio_genes))
        next
      res_prio <- genesAffected[, colnames(genesAffected) %in% prio_genes$Gene]
      res_prio_burden <- as.data.frame(apply(res_prio,1,sum))
      res_prio_burden$SAMPLEID <- genesAffected$SAMPLEID
      anomalyColumn <- paste0(anomaly,"_",string_normalize(term), sep="")
      anomalyColumn <- string (anomalyColumn)

      colnames(res_prio_burden) <- c(anomalyColumn,"SAMPLEID")

      sampleSheetResult <- read.csv2(file.path(resultFolderData,"sample_sheet_result.csv"))
      sampleSheetResult <- sampleSheetResult[, !(colnames(sampleSheetResult)%in% anomalyColumn)]
      sampleSheetResult <- merge(sampleSheetResult, res_prio_burden, by.x="Sample_ID",by.y="SAMPLEID")
      write.csv2(sampleSheetResult,file.path(resultFolderData,"sample_sheet_result.csv"),row.names = FALSE)
    }
  }
}


dir_check_and_create <- function  (baseFolder, subFolders)
{

  parts <- unlist(strsplit(baseFolder, .Platform$file.sep))
  # do.call(file.path, as.list(c(parts[1:length(parts) - 1], subFolders)))
  subFolders <- as.vector(sapply(subFolders, as.character))
  subFolders <-c(parts[1:length(parts)], subFolders)

  for(subFolder in subFolders)
  {
    browser
    if(!exists("subF"))
      subF <- file.path(subFolder)
    else
   {
     subF <- file.path(subF, subFolder)
     if(!dir.exists(subF))
      dir.create(subF)
    }
  }
  return(subF)
}

file_path_build <- function(baseFolder, detailsFilename, extension){

  detailsFilename <- as.vector(sapply(detailsFilename, as.character))
  fileName <- paste0( paste0(detailsFilename, collapse = "_", sep=""),".",extension, sep="")
  fileName <- gsub("__","_", fileName)

  file.path(baseFolder, fileName)

}


string_normalize <- function (string)
{
  string <- gsub("__","_", string)
  string <- gsub(" ", "_", string)
  string <- (gsub("-", "_", string))
  string <- (gsub(":", "_", string))
  string <- (gsub("/", "_", string))
  string <- (gsub("'", "_", string))
  return(toupper(string))

}
