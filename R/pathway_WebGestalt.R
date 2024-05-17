pathway_WebGestalt <- function(study,
  types=c("BP","MF"),  enrich_methods = c("ORA"),
  adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_details,significance,sql_condition="")
{

  #
  # start_fresh <- FALSE
  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  ssEnv <- get_session_info()

  keys <- unique(ssEnv$keys_for_pathway)
  path <- dir_check_and_create(ssEnv$result_folderPathway,"WebGestalt")
  tmp <- tempdir()
  tempFolder <- dir_check_and_create(tmp,c("/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")))


  #check if optional package is installed
  if(!requireNamespace("WebGestaltR", quietly = TRUE))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " WebGestaltR package is not installed. Please install WebGestaltR package to use this function")
    return()
  }

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:(nrow(keys)*length(types)*length(enrich_methods)))
  else
    progress_bar <- ""

  for (i in 1:nrow(keys))
  {
    for( t in 1:length(types))
    {
      type <- types[t]
      for ( em in 1:length(enrich_methods))
      {
        projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=""  ,
          pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)
        filenameResult = file_path_build(path,projectName,"csv")
        if(file.exists(filenameResult))
          next

        #

        enrich_method <- enrich_methods[em]
        if(ssEnv$showprogress)
          progress_bar(sprintf("Searching for disease using webgestalt: %s with %s and %s",keys[i,]$COMBINED,enrich_method,type))
        key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")

        results_inference <- get_results_areas_inference(
          inference_details =  inference_details,
          marker = keys[i,"MARKER"],
          adjust_per_area= adjust_per_area,
          adjust_globally = adjust_globally,
          pvalue_column=  pvalue_column,
          adjustment_method= adjustment_method,
          sql_condition = sql_condition)

        #

        if (nrow(results_inference)==0)
          next

        if(keys[i,]$SUBAREA=="ALL_SUBAREAS")
          gene_set <- results_inference[results_inference$SUBAREA!="WHOLE",]
        else
          gene_set <- results_inference[results_inference$SUBAREA==keys[i,]$SUBAREA,]

        if(keys[i,]$FIGURE=="HYPER_HYPO")
          gene_set <- results_inference[results_inference$FIGURE=="HYPER" | results_inference$FIGURE=="HYPO",]
        else
          gene_set <- results_inference[results_inference$FIGURE==keys[i,]$FIGURE,]

        if(nrow(gene_set)==0)
        {
          log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " No genes found for the key: ", keys[i,])
          next
        }
        # remove duplicates
        gene_set <- aggregate(gene_set[,pvalue_column], by = list(gene_set$AREA_OF_TEST), max)

        colnames(gene_set) <- c("AREA_OF_TEST",pvalue_column)
        gene_set <- gene_set[!duplicated(gene_set$AREA_OF_TEST),]
        gene_set <- gene_set[order(gene_set[,pvalue_column]),]
        gene_set <- na.omit(gene_set)

        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])
        if (nrow(gene_set)==0)
          next
        projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix= paste(type,enrich_method, sep="_") , pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)

        # entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.vector(gene_set$AREA_OF_TEST),column = "ENTREZID", keytype = "SYMBOL")
        # entrez <- unique(entrez)
        # entrez <- na.omit(entrez)

        # log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set (ENTREZ): ",nrow(gene_set), " key: ", keys[i,])
        if (length(gene_set$AREA_OF_TEST)==0)
          next

        geneFile <- file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt")
        write.table(unique(gene_set$AREA_OF_TEST),geneFile,row.names=FALSE,col.names = FALSE,quote =FALSE)

        enrichDataBase <- switch(
          type,
          "BP"="geneontology_Biological_Process",
          "MF"="geneontology_Molecular_Function"
        )

        tryCatch({
          enrichResult <- WebGestaltR::WebGestaltR(
            enrichMethod=enrich_method,
            organism="hsapiens",
            enrichDatabase= enrichDataBase,
            interestGeneFile=geneFile,
            interestGeneType="genesymbol",
            referenceSet= "genome",
            isOutput=FALSE,
            outputDirectory= tempFolder,
            projectName= projectName,
            sigMethod ="top",
            fdrMethod = "BH",
            fdrThr = 1,
            topThr = 200,
            maxNum = 2000,
            minNum = 5
          )
        },
          catch = function(e) {
            log_event(paste("Error in WebGestaltR: ",e))
            return()
          }
        )

        enrichResult$alpha <- as.numeric(ssEnv$alpha)
        enrichResult$pvalue_column <- pvalue_column
        enrichResult$type <- type
        enrichResult$enrich_method <- enrich_method
        enrichResult$key <- key
        enrichResult$MARKER <- keys[i,]$MARKER
        enrichResult$FIGURE <- keys[i,]$FIGURE
        enrichResult$AREA <- keys[i,]$AREA
        enrichResult$SUBAREA <- keys[i,]$SUBAREA

        # format FDR to use scientific notation

        enrichResultToPlot <- enrichResult[1:15,]
        # browser
        # if(min(enrichResultToPlot$FDR) < as.numeric((ssEnv$alpha)))
        # {
        plotFileName = file_path_build(path,projectName,ssEnv$plot_format)
        grDevices::png(file= plotFileName, width=2048,height=2048, bg = "transparent")
        enrichResultToPlot <- as.data.frame(enrichResultToPlot)
        enrichResultToPlot <- enrichResultToPlot[order(enrichResultToPlot$expect),]
        graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
        par(cex.lab = 4)
        # Create a color vector based on the condition FDR < 0.05
        bar_colors <- ifelse(enrichResultToPlot$FDR < as.numeric((ssEnv$alpha)), ssEnv$color_palette[1], ssEnv$color_palette[2])
        # Generate the barplot with colors
        graphics::barplot(height = enrichResultToPlot$expect,
          names = enrichResultToPlot$description,
          # names = wrap_it(paste0(enrichResultToPlot$description,
          #   " (FDR=", format(enrichResultToPlot$FDR, scientific = TRUE), ")",
          #   " (Enr. Ratio = ", round(enrichResultToPlot$enrichmentRatio, 2), ")",
          #   sep=" "), 25),
          horiz = TRUE,
          las = 1,
          cex.names = 2,
          space = 0.1,
          col = bar_colors) # Use the color vector here
        # graphics::barplot(height=enrichResultToPlot$expect, names=wrap_it(paste0(enrichResultToPlot$description,"(FDR=", round(enrichResultToPlot$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
        graphics::mtext( paste0(keys[i,"AREA"]," ",keys[i,"SUBAREA"]," ", keys[i,"MARKER"], " ",keys[i,"FIGURE"]," ",pvalue_column,sep=" "), side = 3, line = 1, cex = 2)
        grDevices::dev.off()
        # }
        file.remove(geneFile)

        unlink(geneFile)

        if(exists("enrichResult"))
          if(exists("enrichResultFinal"))
            enrichResultFinal <- rbind(enrichResultFinal,enrichResult)
        else
          enrichResultFinal <- enrichResult

      }
    }

    if(exists("enrichResultFinal"))
    {
      projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=""  , pvalue_column=pvalue_column, as.numeric(ssEnv$alpha), significance)
      filenameResult = file_path_build(path,projectName,"csv")
      write.csv2(enrichResultFinal, filenameResult)
      rm(enrichResultFinal)
    }
  }
}


# Core wrapping function
wrap_it <- function(x, len)
{
  sapply(x, function(y) paste0(strwrap(y, len),
    collapse = "\n"),
    USE.NAMES = FALSE)
}



