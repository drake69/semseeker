pathway_WebGestalt <- function(study,
  types=c("BP","MF"),  enrich_methods = c("ORA"),
  pvalue = 0.05, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_details,result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{

  start_fresh <- FALSE
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)

  keys <- unique(ssEnv$keys_for_pathway)
  path <- semseeker:::dir_check_and_create(ssEnv$result_folderPathway,"WebGestalt")
  tmp <- tempdir()
  tempFolder <- semseeker:::dir_check_and_create(tmp,c("/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")))


  #check if optional package is installed
  if(!requireNamespace("WebGestaltR", quietly = TRUE))
  {
    log_event("WebGestaltR package is not installed. Please install pathfindR package to use this function")
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
        projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=""  , pvalue_column=pvalue_column, pvalue)
        filenameResult = semseeker:::file_path_build(path,projectName,"csv")
        if(file.exists(filenameResult))
          next


        enrich_method <- enrich_methods[em]
        if(ssEnv$showprogress)
          progress_bar(sprintf("Searching for disease using webgestalt: %s with %s and %s",keys[i,]$COMBINED,enrich_method,type))
        key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")
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
        projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix= paste(type,enrich_method, sep="_") , pvalue_column=pvalue_column, pvalue)

        geneFile <- file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt")
        write.table(unique(results_inference$AREA_OF_TEST),geneFile,row.names=FALSE,col.names = FALSE,quote =FALSE)

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
            sigMethod ="top"
          )

          enrichResult$pvalue <- pvalue
          enrichResult$pvalue_column <- pvalue_column
          enrichResult$type <- type
          enrichResult$enrich_method <- enrich_method
          enrichResult$key <- key

          if(min(enrichResult$FDR) < as.numeric(ssEnv$alpha))
          {
            plotFileName = semseeker:::file_path_build(path,projectName,"png")
            grDevices::png(file= plotFileName, width=2048,height=2048, bg = "transparent")
            enrichResult <- as.data.frame(enrichResult)
            enrichResult <- enrichResult[order(enrichResult$expect),]
            graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
            par(cex.lab = 4)
            graphics::barplot(height=enrichResult$expect, names=wrap_it(paste0(enrichResult$description,"(FDR=", round(enrichResult$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
            graphics::mtext( paste0(keys[i,"AREA"],keys[i,"SUBAREA"], keys[i,"MARKER"], keys[i,"FIGURE"],pvalue_column,sep=" "), side = 3, line = 1, cex = 2)
            grDevices::dev.off()
          }
          file.remove(geneFile)
        },
          catch = function(e) {
            log_event(paste("Error in WebGestaltR: ",e))
            return()
          }
          )

        if(exists("enrichResult"))
          if(exists("enrichResultFinal"))
            enrichResultFinal <- rbind(enrichResultFinal,enrichResult)
          else
            enrichResultFinal <- enrichResult

      }
    }

    if(exists("enrichResultFinal"))
    {
      projectName <- phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=""  , pvalue_column=pvalue_column, pvalue)
      filenameResult = semseeker:::file_path_build(path,projectName,"csv")
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



