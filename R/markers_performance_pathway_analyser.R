#' @export
markers_performance_pathway_analyser <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",
  significance = TRUE, ...)
{
  # #
  # inference_details <- expand.grid(
  #   "independent_variable"= c("Sample_Group"),
  #   "covariates"=c(""),
  #   "family_test"=c("wilcoxon"),
  #   "transformation"=c("scale","none"),
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE
  # )
  # pvalue_column="PVALUE_ADJ_ALL_BH"
  # result_folder <- "~/Documents/Dati_Lavoro/pediatric_leukemia/results/BALLBM/"
  # ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE)


  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  keys <- unique(ssEnv$keys_for_pathway)
  inference_details <- as.data.frame(inference_details)
  key_pathway <- data.frame("type"="Pathway", "label"="pathfindR","column_of_id"="ID","column_of_description"="Term_Description", "column_of_pvalue"="highest_p","column_of_enrichment"="Fold_Enrichment")
  key_pathway <- rbind(key_pathway, data.frame("type"="Pathway","label"="WebGestalt","column_of_id"="geneSet","column_of_description"="description","column_of_pvalue"="FDR","column_of_enrichment"="enrichmentRatio"))
  key_pathway <- rbind(key_pathway, data.frame("type"="Phenotype","label"="phenolyzer","column_of_id"="Description","column_of_description"="HPO","column_of_pvalue"="Score","column_of_enrichment"="Score"))
  #
  for (pt in 1:nrow(key_pathway))
  {
    if(key_pathway[pt,"type"]=="Pathway")
      path <- dir_check_and_create(ssEnv$result_folderPathway,key_pathway[pt,"label"])
    else
      path <- dir_check_and_create(ssEnv$result_folderPhenotype,key_pathway[pt,"label"])

    column_of_id <- key_pathway[pt,"column_of_id"]
    column_of_pvalue <- key_pathway[pt,"column_of_pvalue"]
    column_of_description <- key_pathway[pt,"column_of_description"]
    column_of_enrichment <- key_pathway[pt,"column_of_enrichment"]

    for (id in 1:nrow(inference_details))
    {
      # id <- 1
      # id <- 2
      inference_detail <- inference_details[id,]
      family_test <- inference_detail$family_test
      transformation <- as.character(inference_detail$transformation)
      file_prfx <- paste(pvalue_column,ssEnv$alpha,transformation, family_test,sep="_")
      if (exists("aggregated_patwhay_result_total"))
        rm(aggregated_patwhay_result_total)

      for (i in 1:nrow(keys))
      {
        # i <- 1
        # i <- 2
        file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",suffix="",
          pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
        file_name = file_path_build(path,file_name,"csv")
        # #
        if (!file.exists(file_name))
        {
          log_event("DEBUG: patwhay_result ", file_name, " is missed !")
          next
        }
        # read the patwhay_result
        patwhay_result <- read.csv2(file_name)
        # patwhay_result$key <- paste(patwhay_result$MARKER,patwhay_result$FIGURE,patwhay_result$AREA,patwhay_result$SUBAREA,sep="_")
        if (!exists("aggregated_patwhay_result_total"))
          aggregated_patwhay_result_total <- patwhay_result
        else
          aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, patwhay_result)
      }

      if (!exists("aggregated_patwhay_result_total"))
        next
      if(nrow(aggregated_patwhay_result_total) == 0)
        next

      aggregated_patwhay_result <- aggregated_patwhay_result_total
      #
      fdr <- aggregate(aggregated_patwhay_result[, column_of_pvalue], by = list(aggregated_patwhay_result[,column_of_id]), FUN = max)
      colnames(fdr) <- c( column_of_id, column_of_pvalue)
      # #
      aggregated_patwhay_result <- unique(aggregated_patwhay_result[,c(column_of_id,"key",column_of_description)])
      aggregated_patwhay_result <- na.omit(aggregated_patwhay_result)
      categories <- unique(na.omit(aggregated_patwhay_result$key))
      # remove _
      categories <- gsub("_", " ", categories)
      if(length(categories)==1)
        next
      split <- split(aggregated_patwhay_result[,column_of_id], aggregated_patwhay_result$key)
      filename <- paste(path, "/",file_prfx,"_venn_diagramm_",key_pathway[pt,"label"],".png",sep = "")

      # overlaps <- Reduce(intersect, split)
      # if(length(overlaps)>0)
      # {
      #   # VennDiagram::venn.diagram(
      #   #   x = split,
      #   #   # fill = ssEnv$color_palette[1:length(split)],
      #   #   fill = rep( x="white", length(split)),
      #   #   alpha = 0.5,
      #   #   category.names = categories,
      #   #   cat.col = rep( x="black", length(split)),
      #   #   cat.cex = 1.2,
      #   #   cat.fontface = "bold",
      #   #   cex = 1.5,
      #   #   output = TRUE,
      #   #   individuals.in.intersections = TRUE,
      #   #   disable.logging = TRUE,
      #   #   filename = filename
      #   # )
      #
      #   # overlaps$key <- aggregated_patwhay_result$key
      #   filename <- paste(path, "/",file_prfx,"_overlaps_webgestalt.csv",sep = "")
      #   write.csv2(overlaps, filename)
      # }


      # aggregated_patwhay_result$key <- gsub("_", " ", aggregated_patwhay_result$key)
      # change colname of column_of_id to column_of_id
      colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result)==column_of_id)] <- "column_of_id_label"
      key_gene_set_pivot <- reshape2::dcast(aggregated_patwhay_result, column_of_id_label ~ key , value.var = "column_of_id_label", fun.aggregate = length)
      key_gene_set_pivot_summary <- reshape2::dcast(aggregated_patwhay_result, column_of_id_label ~ MARKER , value.var = "column_of_id_label", fun.aggregate = length)
      #change back colname of column_of_id to column_of_id
      colnames(key_gene_set_pivot)[which(colnames(key_gene_set_pivot)=="column_of_id_label")] <- column_of_id
      colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result)=="column_of_id_label")] <- column_of_id

      # merge with geneSet name
      tt <- unique(aggregated_patwhay_result[,c(column_of_id,column_of_description)])
      tt <- merge(tt, fdr, by=column_of_id)
      key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary,fdr, by=column_of_id)
      #
      key_gene_set_pivot <- merge(tt,key_gene_set_pivot)
      key_gene_set_pivot_summary <- merge(tt,key_gene_set_pivot_summary)
      # add a column with the total number of gene sets
      key_gene_set_pivot$total <- rowSums(key_gene_set_pivot[,4:ncol(key_gene_set_pivot)])

      #
      # foreach key add a column with column_of_enrichment
      kk <- unique(aggregated_patwhay_result$key)
      for (key in 1:length(kk))
      {
        pp <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$key == kk[key],c(column_of_id,column_of_enrichment)]
        colnames(pp) <- c(column_of_id, paste("enrichment_",kk[key],sep = ""))
        key_gene_set_pivot <- merge(key_gene_set_pivot, pp, by=column_of_id, all.x = TRUE)
      }
      key_gene_set_pivot <- key_gene_set_pivot[order(key_gene_set_pivot$total, decreasing = FALSE),]

      filename <- paste(path, "/",file_prfx,"_pivot_",key_pathway[pt,"label"],".csv",sep = "")
      write.csv2(key_gene_set_pivot, filename)

      filename <- paste(path, "/",file_prfx,"_pivot_summary_",key_pathway[pt,"label"],".csv",sep = "")
      write.csv2(key_gene_set_pivot_summary, filename)

      # key_gene_set_pivot <- subset(key_gene_set_pivot, total < 3)[, colSums(key_gene_set_pivot[,4:ncol(key_gene_set_pivot)]) > 0,]
      # filename <- paste(path, "/",file_prfx,"_reduced_pivot_",key_pathway[pt,"label"],".csv",sep = "")
      # write.csv2(key_gene_set_pivot, filename)

      # plot an heatmap
      # apcluster::heatmap(key_gene_set_pivot[,3:ncol(key_gene_set_pivot)], Rowv=NA, Colv=NA,
      #   col = colorRampPalette(c("white", "blue"))(100), scale="none", margins=c(5,10), main="Gene sets",
      #   xlab="Gene sets", ylab="Key", cexRow=0.5, cexCol=0.5, filename = paste(path, "/",file_prfx,"_heatmap_webgestalt.png",sep = ""))


      # #
      # Apply the function to the split data
      differences <- find_unique_gene_sets(split)
      # find differences
      # differences <- Reduce(setdiff, split)
      diff_df <- data.frame()
      #
      if (length(differences)>0)
      {
        filename <- paste(path, "/",file_prfx,"_differences_",key_pathway[pt,"label"],".csv",sep = "")
        for (d in 1:length(differences))
        {
          tt2 <- merge(differences[d], names(differences[d]) )
          colnames(tt2) <- c(column_of_id,"KEY")
          diff_df <- plyr::rbind.fill(diff_df, tt2)
        }
        diff_df <- merge(diff_df, tt[,c(column_of_id, column_of_description)], by=column_of_id)
        write.csv2(diff_df, filename)
      }
    }

    # remove all files ending with .log
    # files <- list.files(ssEnv$result_folderPathway, pattern = ".log")
    # for ( f in 1:length(files))
    # {
    #   file <- files[f]
    #   unlink(paste(ssEnv$result_folderPathway, "/", file, sep = ""), recursive = TRUE, force = TRUE)
    # }

  }

}

# Function to find unique gene sets for each key
find_unique_gene_sets <- function(split_list) {
  unique_sets <- list()
  keys <- names(split_list)

  for (k in 1:length(keys)) {

    key <- keys[k]
    # Start with the current key's gene sets
    current_sets <- split_list[[key]]

    # Get all other keys
    other_keys <- setdiff(keys, key)

    # Combine gene sets from all other keys
    other_sets <- unlist(split_list[other_keys], use.names = FALSE)

    # Find gene sets unique to the current key
    unique_to_current <- setdiff(current_sets, other_sets)

    # Store the unique gene sets with the key as the name
    unique_sets[[key]] <- unique_to_current
  }

  # remove empty sets
  unique_sets <- unique_sets[sapply(unique_sets, length) > 0]

  return(unique_sets)
}


