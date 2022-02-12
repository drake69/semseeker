createMultipleBed <- function(sampleSheet){

  #create multiple file bed
  # browser()
  keys <- expand.grid("POPULATION"=unique(sampleSheet$Sample_Group),"FIGURE"=keys.figures,"ANOMALY"= keys.anomalies,"EXT"="bed")
  keys <- rbind(keys, expand.grid("POPULATION"=unique(sampleSheet$Sample_Group),"FIGURE"="METHYLATION","ANOMALY"="DELTAS" ,"EXT"="bedgraph"))
  foreach::foreach(i = 1:nrow(keys), cl = computationCluster) %dopar% {
  # for (i in 1:nrow(keys)) {
    # i <- 20
    key <- keys[i,]
    for(j in 1:nrow(sampleSheet))
    {
      # j <- 54
      sample <- sampleSheet[j,]
      # if(sample$Sample_ID=="R07C01_203991450116")
      #    message(j)
      tempresultFolderData <-dir_check_and_create(resultFolderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
      fileToRead <- file_path_build(tempresultFolderData, c(sample$Sample_ID, as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
      message("createMultipleBed file 2 read:", fileToRead)
      if(file.exists(fileToRead))
      {
        message("createMultipleBed read file:", fileToRead)
        localtemp <- utils::read.csv2(fileToRead, sep="\t", header = FALSE)
        localtemp$Sample_ID <- sample$Sample_ID
        if(!exists("localFileRes"))
          localFileRes <- localtemp
        else
          localFileRes <- rbind(localFileRes, localtemp)
      }
    }
    if(exists("localFileRes"))
    {
      fileToWrite <- file_path_build(tempresultFolderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
      utils::write.table(localFileRes, fileToWrite, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      rm(localFileRes)
      message("createMultipleBed, file saved!")
    }
  }
  gc()
}
