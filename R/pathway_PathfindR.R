pathway_pathfindR <- function(study,
  path_db,  iterations = 10, exclude_beta= TRUE,
  pvalue = 0.05, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_details,result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")

  start_fresh <- FALSE
  ssEnv <- init_env(result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  keys <- unique(ssEnv$keys_for_pathway)



  #check if optional package is installed
  if(!requireNamespace("pathfindR", quietly = TRUE))
  {
    log_event("pathfindR package is not installed. Please install pathfindR package to use this function")
    return()
  }

  total_progress <- nrow(keys)*length(path_db)*nrow(inference_details)
  progress <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:total_progress)
  else
    progress_bar <- ""

  for(id in 1:nrow(inference_details))
  {
    inference_detail <- inference_details[id,]
    for(i in 1:nrow(keys))
    {
      seq <- 0
      for(k in 1:length(path_db))
      {
        progress <- progress + 1
        db <- path_db[k]
        # i <- 13
        # db <- "KEGG"
        # pvalue_column <- "PVALUE_ADJ_ALL_BH"

        if(exists("existing_db"))
          existing_keys <- unique(pathway_report[pathway_report$source==db,"key"])

        if(ssEnv$showprogress)
          progress_bar(sprintf("Searching for pathway using pathfindR for %s",keys[i,]$COMBINED))
        results_inference <- get_results_areas_inference(inference_detail,keys[i,]$MARKER, pvalue, adjust_per_area, adjust_globally,pvalue_column,adjustment_method)
        if (nrow(results_inference)==0)
          next


        key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")
        if (nrow(results_inference)==0)
          next
        if(exists("existing_db") & exists("existing_keys"))
          if(any(key %in% existing_keys))
            next
        suffix <- ""
        if(exclude_beta)
          suffix = "without_signal_"

        phenotype_analysis_name <- phenotype_analysis_name(inference_detail, keys[i,],prefix ="", suffix= suffix , pvalue_column, pvalue)
        path <- dir_check_and_create(ssEnv$result_folderPathway,"pathfindR")
        pathway_report_path <- file_path_build(path,phenotype_analysis_name,"csv")

        if(file.exists(pathway_report_path))
        {
          next
          pathway_report <- read.csv2(pathway_report_path)
          existing_db <-   unique(pathway_report$source)
        }

        if(keys[i,]$SUBAREA=="WHOLE")
          gene_set <- results_inference[results_inference$SUBAREA==keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
        else
          gene_set <- results_inference[results_inference$SUBAREA!=keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]

        if(keys[i,]$FIGURE=="BOTH" | keys[i,]$FIGURE=="MEAN")
          gene_set <- results_inference[results_inference$FIGURE==keys[i,]$FIGURE,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
        else
          gene_set <- results_inference[results_inference$FIGURE!="BOTH" & results_inference$FIGURE!="MEAN",c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]

        gene_set <- gene_set[,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column)]
        if(nrow(gene_set)<2)
          next
        if(exclude_beta)
          gene_set$STATISTIC_PARAMETER <- NA

        try(
          {
            output_temp <- pathfindR::run_pathfindR( gene_set , path_db[k], max_gset_size = nrow(gene_set), iterations = iterations,output_dir = paste( tempFolder,"/pathfindR/",study, sep=""),
              plot_enrichment_chart = F)
            seq <- seq + 1

            if(nrow(output_temp)==0)
              next()
            output_temp$key <- key
            output_temp$seq <- seq
            output_temp$gene_count <- nrow(gene_set)
            output_temp$source <- path_db[k]
            output_temp$order <- 1:nrow(output_temp)
          }
        )

        if(exists("result_pathway"))
          result_pathway <- rbind(result_pathway,output_temp)
        else
          result_pathway <- output_temp

      }
      if(exists("result_pathway"))
      {
        if(nrow(result_pathway)!=0)
        {
          result_pathway$PHENOTYPE <- grepl(study,result_pathway$Term_Description, ignore.case = T)
          write.csv2(result_pathway, pathway_report_path)
          rm(result_pathway)
        }
      }
    }

  }
}
