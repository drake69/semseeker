#' @importFrom doRNG %dorng%
create_deltaq <- function(envir, resultPopulation){

  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "envir", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
                           ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                           ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
                           ".RNGkind_length", "tail", "RNGstr","update_multiple_bed")
  i <- 0
  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
                           "FIGURE"=c("BOTH"),
                           "ANOMALY"="DELTAS" ,
                           "EXT"="fst"
  )

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), as.character(key$EXT))
    if(file.exists(fileToRead))
    {
      deltaq_temp <- fst::read.fst(fileToRead, as.data.table = T, columns = c("VALUE"))
    }
    if(exists("deltaq"))
      deltaq <- c(deltaq, deltaq_temp)
    else
      deltaq <- deltaq_temp
  }

  rm(deltaq_temp)
  #value is in the 4th position
  # give to each quantile an even weight
  deltaq <- data.frame("VALUE"=deltaq,"DELTAQ"=as.numeric(dplyr::ntile(x=as.numeric(deltaq) , n=4)) * 2)

  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
                           "FIGURE"=envir$keys_figures_default[,1],
                           "ANOMALY"="DELTAS" ,
                           "EXT"="fst"
  )

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)),  as.character(key$EXT))
    if(file.exists(fileToRead))
    {
      localFileRes <- fst::read.fst(fileToRead, as.data.table = T)
    }

    if(exists("localFileRes"))
      if(!plyr::empty(localFileRes))
      {
        if(ncol(localFileRes)==4)
          colnames(localFileRes) <- c("CHR","START","END","SAMPLEID")
        else
          colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")

        xx <- merge(localFileRes,deltaq, by="VALUE")
        localFileRes[,4] <- xx$DELTAQ

        tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste("DELTAQ_",as.character(key$FIGURE),sep="")))
        fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"), as.character(key$FIGURE)),  as.character(key$EXT))
        update_multiple_bed( fileToWrite, localFileRes)

        tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",
                                         fun.aggregate = sum, drop = TRUE)
        colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
        lbl <- rep(paste("DELTAQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
        data.frame("LABEL"=lbl, tempDataFrame)

        if(exists("deltaq_summary"))
          deltaq_summary <- rbind(deltaq_summary,data.frame("LABEL"=lbl, tempDataFrame))
        else
          deltaq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
      }
  }

  if(!plyr::empty(deltaq_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltaq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
