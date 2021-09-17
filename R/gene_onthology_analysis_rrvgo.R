
geneontology_analysis_rrvgo <- function(resultFolder, fileName){

  localFileName <- file.path(resultFolder, "Inference", fileName)

  geneAnnotatedFile <- read.csv(localFileName, sep=";", dec=",")
  resultFolderrrvgo <- file.path(resultFolder,"Inference","geneontology_rrvgo")
  dir.create(resultFolderrrvgo)
  geneAnnotatedFile <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL" & FIGURE=="BOTH" & (PVALUE  < 0.05 | PVALUEADJ < 0.05))
  # WebGestaltR::idMapping()
  geneAnnotatedFile$AREA_OF_TEST <- biomartr::getGO(genes =  as.vector(geneAnnotatedFile$AREA_OF_TEST), filters  = "ensembl_gene_id", organism="Homo sapiens")
  browser()

  for (pval in c("PVALUE", "PVALUEADJ"))
  {
    geneAnnotated <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL" & FIGURE=="BOTH" & geneAnnotatedFile[,pval]   <0.05)

    keys <- unique(geneAnnotated[, c("ANOMALY","FIGURE","SUBGROUP")] )

    if (nrow(keys)>0)
    for (i in 1:nrow(keys))
    {
        browser()
        key <- keys[i,]
        geneMutation <- subset(geneAnnotated, ANOMALY==key[1,1] & SUBGROUP==key[1,3] & FIGURE==key[1,2])
        tryCatch(
        {
            go_analysis <- subset(geneAnnotated, geneAnnotated$ANOMALY== key[1,1] & geneAnnotated$FIGURE== key[1,2] & geneAnnotated$SUBGROUP == key[1,3] )

            simMatrix <- rrvgo::calculateSimMatrix(go_analysis$AREA_OF_TEST,orgdb="org.Hs.eg.db",ont="BP" ,method="Rel")

            scores <- rrvgo::setNames(go_analysis[, pval], go_analysis$AREA_OF_TEST)

            reducedTerms <- reduceSimMatrix(simMatrix,scores,threshold=0.7,orgdb="org.Hs.eg.db")

            #Scatter plot depicting groups and distance between terms
            localFileName <- paste(folder,"scatterplot_", key[1,1] ,"_",key[1,2] ,"_", key[1,3] ,".png", sep = "")
            ggplot2::scatterPlot(simMatrix, reducedTerms)
            grDevices::dev.off()

            localFileName <- paste(folder,"heatmap_", key[1,1] ,"_",key[1,2] ,"_", key[1,3] ,".png", sep = "")
            # #Similarity matrix heatmap
            grDevices::png(file= localFileName, width=2000, height=2000)
            ggplot2::heatmapPlot(simMatrix,reducedTerms,annotateParent=TRUE,annotationLabel="parentTerm",fontsize=6)
            grDevices::dev.off()


            # #treemap ploty
            localFileName <- paste(folder,"treemap_", key[1,1] ,"_",key[1,2] ,"_", key[1,3] ,".png", sep = "")
            grDevices::png(file= localFileName, width=2000, height=2000)
            ggplot2::treemapPlot(reducedTerms)
            grDevices::dev.off()

        } ,
        finally={
            next
        }
        )
    }
  }
}



