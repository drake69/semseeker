euristic_analysis_webgestalt <- function(resultFolder)
{

   init_env(resultFolder)

  geneAnnotated <- utils::read.csv(file.path(resultFolderData, "/GENE_ANNOTATED.csv"))
  resultFolderDataWebGestalt <- dir_check_and_create(resultFolderDataEuristic ,"Webgestalt")

  geneAnnotated <- subset(geneAnnotated, geneAnnotated$POPULATION != "Reference" & FIGURE =="BOTH")
  geneAnnotated <- unique(geneAnnotated[,c("GENE","GROUP","POPULATION","ANOMALY","FIGURE","SAMPLEID")])
  geneAnnotated$freq <- 1
  keys <- unique(geneAnnotated[,c("GROUP","ANOMALY","FIGURE")])

  if(nrow(keys)>0)
    for(i in 1:nrow(keys))
    {

      figure= keys[i,"FIGURE"]
      group= keys[i,"GROUP"]
      anomaly= keys[i,"ANOMALY"]
      geneMutation <- subset(geneAnnotated, geneAnnotated$FIGURE==figure & geneAnnotatedANOMALY==anomaly & geneAnnotatedGROUP==group)
      if(nrow(geneMutation)==0)
        next
      geneMutation <- reshape2::dcast(data = geneMutation, POPULATION  ~ GENE, value.var = "freq", sum)
      geneMutation <- as.data.frame(t(geneMutation))
      names(geneMutation) <- geneMutation[1,]
      geneMutation <- geneMutation[-1,]
      if(nrow(geneMutation)==0)
        next
      geneMutation$gene <- rownames(geneMutation)
      geneMutation <- subset(geneMutation, geneMutation$Case > 0 & geneMutation$Control ==0)

      utils::write.table(unique(row.names(geneMutation)),
                  file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"),
                  row.names=FALSE,
                  col.names = FALSE ,
                  quote =FALSE)

      projectName <- paste0(anomaly, figure, group, sep="_")
      geneFile <- system.file("extdata", "interestingGenes.txt", package="WebGestaltR")
      tryCatch(
        {
          enrichResult <- WebGestaltR::WebGestaltR(
            enrichMethod="ORA",
            organism="hsapiens",
            enrichDatabase="geneontology_Biological_Process",
            interestGeneFile=geneFile,
            interestGeneType="genesymbol",
            referenceSet= "genome",
            isOutput=FALSE,
            outputDirectory= resultFolderDataWebGestalt,
            projectName= projectName,
            sigMethod ="top")

          if(min(enrichResult$FDR) < 0.05)
             projectName <- paste0("@", projectName)
          filename = file_path_build(resultFolderDataWebGestalt,projectName,"png")
          filenameResult = file_path_build(resultFolderDataWebGestalt,projectName,"csv")
          utils::write.csv(enrichResult, filenameResult)
          grDevices::png(file= filename, width=2480,height=2480)
          enrichResult <- as.data.frame(enrichResult)
          enrichResult <- enrichResult[order(enrichResult$expect),]
          graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
          graphics::barplot(height=enrichResult$expect, names=wrap.it(paste0(enrichResult$description,"(FDR=", round(enrichResult$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
          graphics::mtext( paste0(keys[i,"ANOMALY"], keys[i,"FIGURE"], keys[i,"GROUP"], sep=" "), side = 3, line = 1, cex = 2)
          grDevices::dev.off()

        } ,
        finally={
          next
        }
      )
    }
}
