pathway_burden_CTD <-  function(resultFolder, ctdDataFrameFile ) {

   init_env(resultFolder)
  HYPO <- NULL
  HYPER <- NULL
  POPULATION <- NULL

  ctdDataFrame <- utils::read.csv2(ctdDataFrameFile, sep=",")
  ctdDataFrame <- as.data.frame(ctdDataFrame)
  colnames(ctdDataFrame) <- string_normalize(colnames(ctdDataFrame))

  populations <- c("Control","Case")

  # figures <- c("HYPO", "HYPER", "BOTH")
  # anomalies <- c("MUTATIONS","LESIONS")
  # subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")

  figures <- c("BOTH")
  anomalies <- c("MUTATIONS")
  subGroups <- c("Whole")

  probesPrefix = "PROBES_Gene_"
  mainGroupLabel =  "GENE"
  subGroupLabel="GROUP"

  finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel)

  finalBed$GENE <- string_normalize(finalBed$GENE)

  # browser()
  if (is.null(finalBed))
    return()

  samples <- utils::read.csv(file_path_build(baseFolder =  resultFolderData,detailsFilename = c("sample","sheet","result"),"csv"), sep=";")

  # foreach::foreach(path in ctdDataFrame$PATHWAY) %dopar%

  for (i in 1:nrow(ctdDataFrame))
  {
    path <- ctdDataFrame[i,"PATHWAY.ID"]
    toSplit <- paste(ctdDataFrame[ctdDataFrame$PATHWAY.ID==path, "ANNOTATED.GENES"], collapse = "")
    candidateGenes <- string_normalize(unlist(stringr::str_split(toSplit  ,"\\|")))
    tempDataFrame <- subset(finalBed, finalBed$FIGURE=="BOTH" & finalBed$ANOMALY=="MUTATIONS" & finalBed$GROUP=="Whole")
    tempDataFrame <- subset(tempDataFrame, tempDataFrame$GENE %in% candidateGenes)
    tempDataFrame$path <- string_normalize(path)
    if(nrow(tempDataFrame)==0)
      next
    tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID ~ path, value.var = "freq", sum)
    if(!exists("result"))
      result <- merge(samples, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID",all.x=TRUE)
    else
      result <- merge(result, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID",all.x=TRUE)
  }

  fileName <- file_path_build (dir_check_and_create(resultFolderData,c("Pathway","CTD")), c("mutations","both","burden","pathway"),"csv")
  utils::write.table(result,fileName,row.names = FALSE,col.names = TRUE ,quote = FALSE,sep =";")


}
