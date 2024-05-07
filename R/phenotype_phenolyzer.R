phenotype_phenolyzer <- function(study,
  disease,phenolyzer_folder_bin,minimum_score = 0.5,
  pvalue = 0.05, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_details,result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{

  # check existence of phenolyzer
  if (!file.exists(phenolyzer_folder_bin))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Phenolyzer not found in ",phenolyzer_folder_bin)
    return()
  }

  tmp <- tempdir()

  start_fresh <- FALSE
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)

  diseases <- disease_get(disease =  disease,vocabulary =  "HPO")

  if (is.null(nrow(diseases)))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No disease found for ",disease)
    return()
  }

  keys <- unique(ssEnv$keys_for_pathway)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:(nrow(keys)))
  else
    progress_bar <- ""

  seq <-0
  for(i in 1:nrow(keys))
  {
    random_string <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")
    tempFolder <- dir_check_and_create(tmp,c("/semseeker/",random_string))

    if(ssEnv$showprogress)
      progress_bar(sprintf("Searching for disease using phenolyzer: %s",keys[i,]$COMBINED))
    key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")

    path <- dir_check_and_create(ssEnv$result_folderPhenotype,"phenolyzer")
    phenotype_analysis_name <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=paste(disease,"_gene_scores",sep=""), pvalue_column=pvalue_column, pvalue)
    phenotype_report_path <- file_path_build(path,phenotype_analysis_name,"csv")
    if(file.exists(phenotype_report_path))
      next

    results_inference <- get_results_areas_inference(inference_details,keys[i,"MARKER"], pvalue, adjust_per_area, adjust_globally,pvalue_column,adjustment_method)
    if (nrow(results_inference)==0)
      next
    if(keys[i,]$SUBAREA=="WHOLE")
      gene_set <- results_inference[results_inference$SUBAREA==keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
    else
      gene_set <- results_inference[results_inference$SUBAREA!=keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]

    if(keys[i,]$FIGURE=="BOTH" | keys[i,]$FIGURE=="MEAN")
      gene_set <- results_inference[results_inference$FIGURE==keys[i,]$FIGURE,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
    else
      gene_set <- results_inference[results_inference$FIGURE!="BOTH" & results_inference$FIGURE!="MEAN",c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]

    if (nrow(results_inference)==0)
      next

    phenotype_analysis_name <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=disease, pvalue_column=pvalue_column, pvalue)
    phenotype_report_path <- file_path_build(path,phenotype_analysis_name,"csv")

    gene_set <- results_inference[,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column)]
    gene_set <- na.omit(gene_set)

    if(nrow(gene_set)<2)
      next

    
    file_term <- file.path(tempFolder, paste0("term_",random_string,".txt"))
    file_genes <- file.path(tempFolder, paste0("genes_",random_string,".txt"))
    write.table(unique(gene_set[,"AREA_OF_TEST"]), file_genes, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(diseases, file_term, quote = FALSE, row.names = FALSE, col.names = FALSE)

    try(
      {
        # cpan App::cpanminus
        # biotools:phenolyzer
        # cpanm Module::Name
        # cd /usr/local/lib/
        # git clone https://github.com/WGLab/phenolyzer.git
        # sudo cpan Bio::OntologyIO
        # sudo cpan Graph::Directed

        # perl /usr/local/lib/phenolyzer/disease_annotation.pl -p ./hypertension.txt  -f --gene ./genes_list.txt  -prediction -phenotype -logistic -out hypertension -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15
        pcommand <-
          paste0(
            "cd ",
            tempFolder ,
            " && perl ", phenolyzer_folder_bin ,"disease_annotation.pl " ,
            " -p " ,file_term,
            " -f --gene " ,file_genes,
            " -wordcloud ",
            " -prediction -phenotype -logistic ",
            " -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 ",
            " -nproc ", ssEnv$parallel$nCore ,
            sep = ""
          )

        null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
        message(pcommand)
        sink(null_device)
        system(pcommand)
        sink()

        annotated_gene_file <- file.path(tempFolder,"out.annotated_gene_list")
        if(!file.exists(annotated_gene_file))
        {
          message("No annotated gene file found")
          next
          }
        annotated_genes <-  utils::read.csv2(annotated_gene_file, sep = "\t")


        annotated_genes_temp <- annotated_genes[annotated_genes$Status=="SeedGene",]

        if(nrow(annotated_genes_temp)==0)
          annotated_genes_temp <- annotated_genes[annotated_genes$Score>=minimum_score,]

        if(nrow(annotated_genes_temp)==0)
          next

        seq <- seq + 1
        if(nrow(annotated_genes)==0)
          next()
        annotated_genes$key <- key
        annotated_genes$seq <- seq
        annotated_genes$gene_count <- nrow(gene_set)
        annotated_genes$order <- 1:nrow(annotated_genes)
        annotated_genes$diseases <- paste(diseases$DISEASE, collapse = ";")
        write.csv2(annotated_genes, phenotype_report_path)

        annotated_gene_file <- file.path(tempFolder,"out.annotated_gene_scores")
        phenotype_analysis_name <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=paste(disease,"_gene_scores",sep=""), pvalue_column=pvalue_column, pvalue)
        phenotype_report_path <- file_path_build(path,phenotype_analysis_name,"csv")

        file.copy(annotated_gene_file, phenotype_report_path, overwrite = TRUE)

        
        # phenotype_analysis_name <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=paste(disease,"_gene_cloud",sep=""), pvalue_column=pvalue_column, pvalue)
        # cloud_gene_file <- file.path(tempFolder,"out.annotated_gene_scores")
        # file.copy(cloud_gene_file, phenotype_report_path, overwrite = TRUE)

        unlink(tempFolder, recursive = TRUE)
      }
    )
  }
}
