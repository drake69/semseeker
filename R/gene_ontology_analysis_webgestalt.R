geneontology_analysis_webgestalt <- function(resultFolder, fileName){


  resultFolderDataWebGestalt <- dir_check_and_create(resultFolderData, c("Inference","geneontology_webgestalt", fileName()))

  fileName <- file.path(resultFolderData, "Inference", fileName)
  geneAnnotatedFile <- utils::read.csv(fileName, sep=";", dec=",")
  geneAnnotatedFile <- subset(geneAnnotatedFile, geneAnnotatedFile$GROUP=="GENE" &
                                geneAnnotatedFile$SUBGROUP!="SAMPLE" & geneAnnotatedFile$AREA_OF_TEST !="TOTAL" &
                                geneAnnotatedFile$SUBGROUP!="TOTAL" & geneAnnotatedFile$FIGURE=="BOTH" &
                                (geneAnnotatedFile$PVALUE  <0.05 | geneAnnotatedFile$PVALUEADJ < 0.05))

  for(type in c("BP","MF"))
  {
    for (pval in c("PVALUE", "PVALUEADJ"))
    {
      geneAnnotated <- subset(geneAnnotatedFile, geneAnnotatedFile[,pval]  <0.05)
      keys <- unique(geneAnnotated[, c("ANOMALY","FIGURE","SUBGROUP")] )

      if (nrow(keys)>0)

        for (i in 1:nrow(keys))
        {
          anomaly <- keys[i,"ANOMALY"]
          subgroup <- keys[i,"SUBGROUP"]
          figure <- keys[i,"FIGURE"]

          geneMutation <- subset(geneAnnotated, geneAnnotated$ANOMALY==anomaly & geneAnnotated$SUBGROUP==subgroup & geneAnnotated$FIGURE==figure )
          print(nrow(geneMutation))

          # print(file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"))
          utils::write.table(unique(geneMutation$AREA_OF_TEST),
                      file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"),
                      row.names=FALSE,
                      col.names = FALSE ,
                      quote =FALSE)

          projectName <- paste0(pval, "gene", anomaly, figure, subgroup,type, sep="_")
          tempGeneFileName <- stringi::stri_rand_strings(1, 16, pattern = "[A-Za-z]")
          geneFile <- system.file("extdata", tempGeneFileName, package="WebGestaltR")


              enrichDataBase <- switch(
                type,
                "BP"="geneontology_Biological_Process",
                "MF"="geneontology_Molecular_Function"
              )
              message(enrichDataBase)
              enrichResult <- WebGestaltR::WebGestaltR(
                enrichMethod="ORA",
                organism="hsapiens",
                enrichDatabase= enrichDataBase,
                interestGeneFile=geneFile,
                interestGeneType="genesymbol",
                referenceSet= "genome",
                isOutput=FALSE,
                outputDirectory= resultFolderDataWebGestalt,
                projectName= projectName,
                sigMethod ="top"
              )

              if(min(enrichResult$FDR) < 0.05)
                projectName <- paste0("@", projectName, sep="")

              filename = file_path_build(resultFolderDataWebGestalt,projectName,"png")
              filenameResult = file_path_build(resultFolderDataWebGestalt,projectName,"csv")
              utils::write.csv2(enrichResult, filenameResult)

              grDevices::png(file= filename, width=2048,height=2048, bg = "transparent")
              enrichResult <- as.data.frame(enrichResult)
              enrichResult <- enrichResult[order(enrichResult$expect),]
              graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
              graphics::par(cex.lab = 4)
              graphics::barplot(height=enrichResult$expect, names=wrap.it(paste0(enrichResult$description,"(FDR=", round(enrichResult$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
              graphics::mtext( paste0(pval, "gene", anomaly, figure, subgroup, sep=" "), side = 3, line = 1, cex = 2)
              grDevices::dev.off()
              # break

          file.remove(geneFile)
        }
    }
  }
}
