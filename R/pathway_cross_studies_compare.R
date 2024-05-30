pathway_cross_studies_compare <- function(result_folder, ...){

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  # library(VennDiagram)
  pathway_inference <- pathway_inference[pathway_inference$occurrence >= occurrence_threshold,]
  studies_to_comb <- na.omit(unique(pathway_inference$STUDY))
  for( j in length(studies_to_comb):1)
  {
    studies_comb <- combinat::combn(studies_to_comb, j)
    if (j==length(studies_to_comb))
      studies_comb <- data.frame("st"=studies_comb)
    log_event(studies_comb)
    for(k in 1: ncol(studies_comb))
    {
      pathway_inference_comb <- subset(pathway_inference, STUDY %in% studies_comb[,k])
      keys <- unique(pathway_inference_comb[, c("key")])
      dest_folder <- paste(result_folder, subFolder, sep="")
      study_count <- length(unique(pathway_inference_comb$STUDY))
      if(!dir.exists(dest_folder))
        dir.create(dest_folder)
      # Load library
      for (i in 1:length(keys))
      {
        pathway_set <- pathway_inference_comb[pathway_inference_comb$key == keys[i],]
        pathway_set <- pathway_set[ ,c("ID","STUDY") ]
        pathway_set <- na.omit(pathway_set)
        categories <- unique(pathway_set$STUDY)
        if(length(categories)==1)
          next
        SPLIT <- split(pathway_set$ID, pathway_set$STUDY)
        studies <- paste(sort(categories), collapse = "-")
        resultFolder <- paste(result_folder, subFolder, studies, sep="/" )
        if(!file.exists(resultFolder))
          dir.create(resultFolder)
        overlaps <- Reduce(intersect, SPLIT)
        if(length(overlaps)>0)
        {
          if(!dir.exists(resultFolder))
            dir.create(resultFolder)
          filename <-
            paste(
              resultFolder,  "/PATHWAY_",
              keys[i],
              "_venn_diagramm.png",
              sep = ""
            )

          # Chart
          # Set up the Venn diagram parameters
          color_palette <- color_palette[length(SPLIT)]
          # Plot the Venn diagram
          VennDiagram::venn.diagram(
            x = SPLIT,
            fill = color_palette,
            alpha = 0.5,
            category.names = categories,
            cat.col = color_palette,
            cat.cex = 1.2,
            cat.fontface = "bold",
            cex = 1.5,
            # main = "Venn Diagram",
            # sub = "Example Venn diagram with custom colors",
            filename = filename
          )

          filename <-
            paste(
              resultFolder,  "/PATHWAY_",
              keys[i],
              "_pathways.csv",
              sep = ""
            )
          write.csv2(overlaps, filename)
          pathway_shared <- pathway_inference_comb[pathway_inference_comb$key == keys[i] & pathway_inference_comb$ID %in% overlaps, ]
          filename <-
            paste(
              resultFolder,  "/PATHWAY_",
              keys[i],
              "_pathways_shared.csv",
              sep = ""
            )
          write.csv2(pathway_shared, filename)
        }
      }
    }
  }

  for (i in 1:length(keys))
  {
    pathway_set <- pathway_inference_comb[pathway_inference_comb$key == keys[i],]
    pathway_set <- pathway_set[ ,c("ID","STUDY") ]
    pathway_set <- na.omit(pathway_set)
    categories <- unique(pathway_set$STUDY)
    if(length(categories)==1)
      next
    SPLIT <- split(pathway_set$ID, pathway_set$STUDY)
    resultFolder <- paste("~/Desktop/VENN", subFolder, sep="/" )
    overlaps <- Reduce(intersect, SPLIT)
    if(length(overlaps)>0)
    {
      if(!dir.exists(resultFolder))
        dir.create(resultFolder)
      filename <-
        paste(
          resultFolder,  "/ALL_PATHWAYS_",
          keys[i],
          "_venn_diagramm.png",
          sep = ""
        )
      # Chart
      VennDiagram::venn.diagram(
        x = SPLIT,
        category.names = categories,
        filename = filename,
        output = TRUE,
        individuals.in.intersections = TRUE,
        disable.logging = T
      )
    }
  }
}
