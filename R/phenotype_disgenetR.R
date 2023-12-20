phenotype_disgenetR <- function(inference_details,marker, study,
  disease , score_limit = c(0.2,1),
  pvalue = 0.05, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUEADJ_ALL_BH")
{

  ssEnv <- get_session_info()
  resultFolder <- ssEnv$result_folderInference
  results_inference <- get_results_inference_genes(inference_details,marker, pvalue, adjust_per_area, adjust_globally,pvalue_column,adjustment_method)

  # require(disgenet2r)
  phenotype_report_path <-  paste(resultFolder,"/Phenotype/",sep="")
  if(!dir.exists(paste(phenotype_report_path)))
    dir.create(phenotype_report_path)


  keys <- unique(results_inference[ ,c("SUBAREA","AREA","MARKER","FIGURE")])

  for(i in 1:nrow(keys))
  {
    key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")
    phenotype_report_path <-  paste(resultFolder,"/Disease/disgenet_", keys[i,]$MARKER,"_",gsub("[.]","",pvalue),"_",transformation,"_",FAMILY,"_",applied_model,"_", adjustment_text,"_phenotype_result.csv", sep="")
    gene_set <- results_inference[
      results_inference$SUBAREA==keys[i,]$SUBAREA &
        results_inference$AREA==keys[i,]$AREA &
        results_inference$MARKER==keys[i,]$MARKER &
        results_inference$FIGURE==keys[i,]$FIGURE
      ,c("AREA_OF_TEST","BETA",pvalue_column,"PVALUE"),]

    gene_set <- subset(gene_set, gene_set[,pvalue_column] < pvalue)
    gene_set <- gene_set[,c("AREA_OF_TEST","BETA",pvalue_column)]
    gene_set <- na.omit(gene_set)
    if(nrow(gene_set)<2)
      next
    try(
      {
        result <- disease2evidence(
          disease  = disease,
          gene = as.vector(gene_set$AREA_OF_TEST),
          database = "ALL",
          score    = score_limit
        )
        output_temp <- extract(result)

        seq <- seq + 1
        if(nrow(output_temp)==0)
          next()
        output_temp$key <- key
        output_temp$seq <- seq
        output_temp$gene_count <- nrow(gene_set)
        output_temp$source <- path_db[k]
        output_temp$order <- 1:nrow(output_temp)
        if(exists("phenotype_result"))
          phenotype_result <- rbind(phenotype_result, output_temp)
        else
          phenotype_result <- output_temp
        write.csv(phenotype_result, phenotype_report_path)
      }
    )
  }
}
