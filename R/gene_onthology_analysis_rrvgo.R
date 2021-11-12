
geneontology_analysis_rrvgo <- function(resultFolder, fileName){

  resultFolderData <- file.path(resultFolder,"Data")
  resultFolderrrvgo <- dir_check_and_create(resultFolderData,"Inference","Geneontology","rrvgo",fileName)
  dataFileName <- file.path(resultFolderData, "Inference", fileName)

  geneAnnotatedFile <- read.csv(dataFileName, sep=";", dec=",")
  geneAnnotatedFile <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL")
  geneAnnotatedFile <- subset(geneAnnotatedFile, FIGURE=="BOTH")
  geneAnnotatedFile <- subset(geneAnnotatedFile,PVALUE  < 0.05 | PVALUEADJ < 0.05)

  subgroups <- unique(geneAnnotatedFile$SUBGROUP)

  # library("GSEABase")
  frame <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO)
  goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  goFrame <- AnnotationDbi::GOFrame(goframeData,organism="Homo sapiens")
  goAllFrame <- AnnotationDbi::GOAllFrame(goFrame)
  gsc <- GSEABase::GeneSetCollection(goAllFrame, setType = GSEABase::GOCollection())
  universe = AnnotationDbi::Lkeys(org.Hs.eg.db::org.Hs.egGO)
  # dbDisconnect()

  for(type in c("BP","MF"))
  {
    # type <- "BP"
    # message(type)
    semdata <- GOSemSim::godata("org.Hs.eg.db", keytype = "ENTREZID", ont=type)
    for (subgroup in subgroups)
    {
      for (pval in c("PVALUE", "PVALUEADJ"))
      {
        # pval <- "PVALUE"
        # message(pval)
        keys <- unique(geneAnnotatedFile[, c("ANOMALY","FIGURE","SUBGROUP")] )

        if (nrow(keys)>0)
          for (i in 1:nrow(keys))
          {
            # i <- 1
            # message(i)
            key <- keys[i,]
            anomaly <- key[1,1]
            figure <- key[1,2]
            subgroup <- key[1,3]

            message(type, " ",pval, " ", anomaly, " ", figure, " ", subgroup, " ")
            geneAnnotated <- subset(geneAnnotatedFile, ANOMALY==anomaly & FIGURE==figure & SUBGROUP== subgroup & geneAnnotatedFile[,pval]   <0.05)
            if(nrow(geneAnnotated) ==0)
              next

            searchedGenes <- gsub("_","-",geneAnnotated$AREA_OF_TEST)
            resultTemp <- AnnotationDbi::mget(x=searchedGenes,envir=org.Hs.egALIAS2EG,mode = "any",ifnotfound = NA,inherits = FALSE)
            foundGenes <- names(resultTemp)
            notfoundGenes <- gsub("_","-",searchedGenes[!( searchedGenes %in% foundGenes )])
            if(length(notfoundGenes)>0)
              browser()
            # resultTemp1 <- AnnotationDbi::mget(x=notfoundGenes,envir=org.Hs.egALIAS2EG,mode = "any",ifnotfound = NA,inherits = FALSE)
            # foundGenes1 <- names(resultTemp1)
            # notfoundGenes1 <- notfoundGenes[!( notfoundGenes %in% foundGenes1 )]

            genes <- unique(as.character(unlist(resultTemp)))

            # genes <- (as.character(unlist(mget(x=gsub("_","-",geneAnnotated$AREA_OF_TEST),envir=org.Hs.egSYMBOL2EG))))

            params <- Category::GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                                   geneSetCollection=gsc,
                                                   geneIds = genes,
                                                   universeGeneIds = universe,
                                                   ontology = type,
                                                   pvalueCutoff = 0.05,
                                                   conditional = FALSE,
                                                   testDirection = "over")

            Over <- GOstats::hyperGTest(params)
            go_analysis <- GOstats::summary(Over)

            # browser()

            if(length(go_analysis)==0 | nrow(go_analysis) < 2)
            {
              next
            }

            colName <- switch(
              type,
              "BP"="GOBPID",
              "MF"="GOMFID"
            )

                simMatrix <- rrvgo::calculateSimMatrix(go_analysis[,colName]  ,orgdb="org.Hs.eg.db", ont=type ,method="Rel", keytype = "ENTREZID", semdata = semdata )
                scores <- setNames(-log10(go_analysis$Pvalue), go_analysis[,colName] )

                reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,scores,threshold=0.7,orgdb="org.Hs.eg.db")


                typeTitle = switch(
                  type,
                  "BP" = "biological process",
                  "MF" = "molecular function"
                )

                figureTitle <- "hypermethylated and hypomethylated"

                #Scatter plot depicting groups and distance between terms
                plotFileName <- paste0(pval,"_",type, "_","scatterplot_", anomaly ,"_",figure ,"_", subgroup ,".png", sep = "")
                # if(min(enrichResult$FDR)<0.05)
                #   plotFileName <- paste0("@", plotFileName, sep="")
                plotFileName <- paste0(resultFolderrrvgo,"/", plotFileName, sep="")
                sp <- rrvgo::scatterPlot(simMatrix, reducedTerms)
                sp <- sp + ggplot2::ggtitle(paste0(typeTitle," ", anomaly," ", figureTitle, " ", subgroup ," ", pval, sep = "")) +ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
                ggplot2::ggsave(plot = sp, filename = plotFileName, width = 10, height = 10, device='png', dpi=600, units = "in", bg="transparent")
                # grDevices::dev.off()


                plotFileName <- paste0(pval,"_",type, "_","heatmap_", anomaly ,"_",figure ,"_", subgroup ," ", pval,".png", sep = "")
                # if(min(enrichResult$FDR)<0.05)
                #   plotFileName <- paste0("@", plotFileName, sep="")
                plotFileName <- paste0(resultFolderrrvgo,"/", plotFileName, sep="")
                # #Similarity matrix heatmap
                grDevices::png(file= plotFileName, width=10, height=10, res = 600, units = "in", bg="transparent")
                rrvgo::heatmapPlot(simMatrix,reducedTerms,annotateParent=TRUE,annotationLabel="parentTerm",fontsize=6)
                # is not possible to add a title !
                # hp <- hp + ggplot2::ggtitle(paste0(typeTitle, " ", anomaly ," ",figureTitle ," ", subgroup , sep = "")) + theme(plot.title = element_text(hjust = 0.5))
                # ggplot2::ggsave(plot = hp, filename = plotFileName, width = 10, height = 10, device='png', dpi=600, units = "in", bg="transparent")
                grDevices::dev.off()

                # #treemap ploty
                # the returned value of treemaplot is not a ggplot2 object is a list!
                plotFileName <- paste0(pval,"_",type, "_","treemap_", anomaly ,"_",figure ,"_", subgroup ," ", pval,".png", sep = "")
                # if(min(enrichResult$FDR)<0.05)
                #   plotFileName <- paste0("@", plotFileName, sep="")
                plotFileName <- paste0(resultFolderrrvgo,"/", plotFileName, sep="")
                grDevices::png(file= plotFileName, width=10, height=10 , res = 600, units = "in", bg="transparent")
                rrvgo::treemapPlot(reducedTerms, title = paste0(typeTitle, " ", anomaly ," ",figureTitle ," ", subgroup ," ", pval, sep = ""))
                # tp <- tp + ggplot2::ggtitle(paste0(typeTitle, " ", anomaly ," ",figureTitle ," ", subgroup , sep = ""))
                # ggsave(plot = tp, filename = plotFileName, width = 10, height = 10, device='png', dpi=600, units = "in", bg="transparent")
                grDevices::dev.off()

          }
      }
    }
  }
}



