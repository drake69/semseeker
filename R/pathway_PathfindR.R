pathway_pathfindR <- function(inference_details,marker,study,
  path_db,  iterations = 10, exclude_beta= FALSE,
  pvalue = 0.05,  adjust_per_area = F, adjust_globally = F,  adjustment_method = "BH", pvalue_column="PVALUEADJ_ALL_BH")
{

  results_inference <- get_results_inference_genes(inference_details,marker, pvalue, adjust_per_area, adjust_globally,pvalue_column,adjustment_method)
  pathway_report_path <-  paste(resultFolder,"/Pathway/",sep="")

  if(!dir.exists(paste(pathway_report_path)))
    dir.create(pathway_report_path)

  for(k in 1:length(path_db))
  {
    db <- path_db[k]
    if(exists("existing_db"))
      existing_keys <- unique(pathway_report[pathway_report$source==db,"key"])
    for(i in 1:nrow(keys))
    {
      # k <- 1
      # i <- 10
      key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")
      if(exists("existing_db") & exists("existing_keys"))
        if(any(key %in% existing_keys))
          next
      if(exclude_beta)
        pathway_report_path <-  paste(resultFolder,"/Pathway/pathfindR_without_beta_", keys[i,]$MARKER,"_",gsub("[.]","",pvalue),"_",transformation,"_",FAMILY,"_",applied_model,"_", adjustment_text,"_pathway_result.csv", sep="")
      else
        pathway_report_path <-  paste(resultFolder,"/Pathway/pathfindR_", keys[i,]$MARKER,"_",gsub("[.]","",pvalue),"_",transformation,"_",FAMILY,"_",applied_model,"_", adjustment_text,"_pathway_result.csv", sep="")
      if(file.exists(pathway_report_path))
      {
        pathway_report <- read.csv2(pathway_report_path)
        existing_db <-   unique(pathway_report$source)
      }
      gene_set <- results_inference[
        results_inference$SUBAREA==keys[i,]$SUBAREA &
          results_inference$AREA==keys[i,]$AREA &
          results_inference$MARKER==keys[i,]$MARKER &
          results_inference$FIGURE==keys[i,]$FIGURE
        ,c("AREA_OF_TEST","BETA",pvalue_column,"PVALUE"),]

      gene_set <- gene_set[,c("AREA_OF_TEST","BETA",pvalue_column)]
      if(nrow(gene_set)<2)
        next
      if(exclude_beta)
        gene_set$BETA <- NA
      try(
        {
          output_temp <- pathfindR::run_pathfindR( gene_set , path_db[k], max_gset_size = nrow(gene_set), iterations = iterations,output_dir = paste("~/Desktop/VENN/report/",study,"/pathfindR", sep=""),
            plot_enrichment_chart = F)
          seq <- seq + 1
          if(nrow(output_temp)==0)
            next()
          output_temp$key <- key
          output_temp$seq <- seq
          output_temp$gene_count <- nrow(gene_set)
          output_temp$source <- path_db[k]
          output_temp$order <- 1:nrow(output_temp)
          if(exists("pathway_result"))
            pathway_result <- rbind(pathway_result, output_temp)
          else
            pathway_result <- output_temp
          write.csv(pathway_result, pathway_report_path)
        }
      )
    }
  }
}

