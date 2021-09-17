geneontology_analysis_webgestalt <- function(resultFolder, fileName){

  fileName <- file.path(resultFolder, "Inference", fileName)

  geneAnnotated <- read.csv(fileName, sep=";", dec=",")
  resultFolderWebGestalt <- file.path(resultFolder,"Inference","geneontology_webgestalt")
  dir.create(resultFolderWebGestalt)
  geneAnnotatedFile <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL" & FIGURE=="BOTH" & (PVALUE  <0.05 | PVALUEADJ < 0.05))

  for (pval in c("PVALUE", "PVALUEADJ"))
  {
    geneAnnotated <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL" & FIGURE=="BOTH" & geneAnnotatedFile[,pval]  <0.05)

    keys <- unique(geneAnnotated[, c("ANOMALY","FIGURE","SUBGROUP")] )

    if (nrow(keys)>0)
      for (i in 1:nrow(keys))
      {
        geneMutation <- subset(geneAnnotated, ANOMALY==keys[i,"ANOMALY"] & SUBGROUP==keys[i,"SUBGROUP"] & FIGURE==keys[i,"FIGURE"] )
        print(nrow(geneMutation))

        print(file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"))
        write.table(unique(geneMutation$AREA_OF_TEST),
                    file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"),
                    row.names=FALSE,
                    col.names = FALSE ,
                    quote =FALSE)

        projectName <- paste(pval, "gene", keys[i,"ANOMALY"], keys[i,"FIGURE"], keys[i,"SUBGROUP"], sep="_")
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
              outputDirectory= resultFolderWebGestalt,
              projectName= projectName,
              sigMethod ="top")

            filename = file.path(paste( resultFolderWebGestalt,"/",projectName,".png",sep=""))
            filenameResult = file.path(paste( resultFolderWebGestalt,"/",projectName,".csv",sep=""))
            write.csv2(enrichResult, filenameResult)
            grDevices::png(file= filename, width=2000, height=2000)
            enrichResult <- as.data.frame(enrichResult)
            enrichResult <- enrichResult[order(enrichResult$expect),]
            graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
            par(cex.lab = 4)
            graphics::barplot(height=enrichResult$expect, names=wrap.it(paste(enrichResult$description,"(FDR=", round(enrichResult$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
            graphics::mtext( paste(pval, "gene", keys[i,"ANOMALY"], keys[i,"FIGURE"], keys[i,"SUBGROUP"], sep=" "), side = 3, line = 1, cex = 2)
            grDevices::dev.off()
            # break
          } ,
          finally={
            # break
            next
          }
        )
        file.remove(geneFile)
      }
  }
}
