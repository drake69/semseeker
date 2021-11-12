geneontology_analysis_webgestalt <- function(resultFolder, fileName){


  resultFolderWebGestalt <- check_and_create_dir(resultFolder, c("Inference","geneontology_webgestalt", fileName()))

  fileName <- file.path(resultFolder, "Inference", fileName)
  geneAnnotatedFile <- read.csv(fileName, sep=";", dec=",")
  geneAnnotatedFile <- subset(geneAnnotatedFile, GROUP=="GENE" & SUBGROUP!="SAMPLE" & AREA_OF_TEST !="TOTAL" & SUBGROUP!="TOTAL" & FIGURE=="BOTH" & (PVALUE  <0.05 | PVALUEADJ < 0.05))

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

          geneMutation <- subset(geneAnnotated, ANOMALY==anomaly & SUBGROUP==subgroup & FIGURE==figure )
          print(nrow(geneMutation))

          # print(file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"))
          write.table(unique(geneMutation$AREA_OF_TEST),
                      file.path(system.file(package="WebGestaltR"),"extdata/interestingGenes.txt"),
                      row.names=FALSE,
                      col.names = FALSE ,
                      quote =FALSE)

          projectName <- paste0(pval, "gene", anomaly, figure, subgroup,type, sep="_")
          tempGeneFileName <- stringi::stri_rand_strings(1, 16, pattern = "[A-Za-z0-9]")
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
                outputDirectory= resultFolderWebGestalt,
                projectName= projectName,
                sigMethod ="top"
              )

              if(min(enrichResult$FDR) < 0.05)
                projectName <- paste0("@", projectName, sep="")

              filename = file_path_build(resultFolderWebGestalt,projectName,"png")
              filenameResult = file_path_build(resultFolderWebGestalt,projectName,"csv")
              write.csv2(enrichResult, filenameResult)

              grDevices::png(file= filename, width=2048,height=2048, bg = "transparent")
              enrichResult <- as.data.frame(enrichResult)
              enrichResult <- enrichResult[order(enrichResult$expect),]
              graphics::par(mar = c(5, 30, 5, 5)) # Set the margin on all sides to 6
              par(cex.lab = 4)
              graphics::barplot(height=enrichResult$expect, names=wrap.it(paste0(enrichResult$description,"(FDR=", round(enrichResult$FDR,2),")", sep=" "),25) , horiz=T, las = 1 ,cex.names= 2, space=0.1)
              graphics::mtext( paste0(pval, "gene", anomaly, figure, subgroup, sep=" "), side = 3, line = 1, cex = 2)
              grDevices::dev.off()
              # break

          file.remove(geneFile)
        }
    }
  }
}
