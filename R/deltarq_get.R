#' @importFrom doRNG %dorng%
deltarq_get <- function( resultPopulation){

  ssEnv <- get_session_info()

  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "ssEnv", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
    ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
    ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
    ".RNGkind_length", "tail", "RNGstr","update_multiple_bed","probes_get")
  i <- 2
  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
    "FIGURE"=ssEnv$keys_figures_default[,1],
    "ANOMALY"="DELTAR" ,
    "EXT"="fst"
  )

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), as.character(key$EXT))
    if(file.exists(fileToRead))
    {
      deltarq_temp <- fst::read.fst(fileToRead, as.data.table = T)
      colnames(deltarq_temp) <- c("CHR","START","END","VALUE","SAMPLEID")
      deltarq_temp$FIGURE <- key$FIGURE
      deltarq_temp$POPULATION <- key$POPULATION
      if(exists("deltarq"))
        deltarq <- rbind(deltarq, deltarq_temp)
      else
        deltarq <- deltarq_temp
    }
  }

  if (!exists("deltarq")| plyr::empty(deltarq))
  {
    stop("Something wrong with multiple bed files!")
  }

  deltarq$DELTARQ <- as.numeric(dplyr::ntile(x=deltarq[,"VALUE"] , n=100))

  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
    "FIGURE"=ssEnv$keys_figures_default[,1],
    "ANOMALY"="DELTAS" ,
    "EXT"="fst"
  )

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]

    if(exists("deltarq"))
    {
      localFileRes <- deltarq[ deltarq$POPULATION==as.character(key$POPULATION)
        & deltarq$FIGURE==as.character(key$FIGURE)
        ,c("CHR","START","END","DELTARQ","SAMPLEID")]
      if(!plyr::empty(localFileRes))
      {
        colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")

        tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$POPULATION) ,paste("DELTARQ_",as.character(key$FIGURE),sep="")))
        fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTARQ"), as.character(key$FIGURE)),  as.character(key$EXT))
        fst::write.fst( x= localFileRes, fileToWrite)
        message("INFO: ", Sys.time(), " Created DELTARQ multiple annotated file!", fileToWrite)

        tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
        colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
        lbl <- rep(paste("DELTARQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
        data.frame("LABEL"=lbl, tempDataFrame)

        if(exists("deltarq_summary"))
          deltarq_summary <- rbind(deltarq_summary,data.frame("LABEL"=lbl, tempDataFrame))
        else
          deltarq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
      }
    }
  }


  # for(i in 1:nrow(localKeys))
  # {
  #   key <- localKeys[i,]
  #   tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
  #   fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)),  as.character(key$EXT))
  #   if(file.exists(fileToRead))
  #   {
  #     localFileRes <- fst::read.fst(fileToRead, as.data.table = T)
  #
  #     if(!plyr::empty(localFileRes))
  #     {
  #       colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")
  #       deltarq <- unique(deltarq)
  #       localFileRes[,4] <- deltarq[ deltarq$VALUE %in% localFileRes$VALUE, "DELTARQ"]
  #
  #       tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$POPULATION) ,paste("DELTARQ_",as.character(key$FIGURE),sep="")))
  #       fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTARQ"), as.character(key$FIGURE)),  as.character(key$EXT))
  #       fst::write.fst( x= localFileRes, fileToWrite)
  #       # localFileRes <-fst::read.fst(fileToWrite)
  #
  #       tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
  #       colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
  #       lbl <- rep(paste("DELTARQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
  #       data.frame("LABEL"=lbl, tempDataFrame)
  #
  #       if(exists("deltarq_summary"))
  #         deltarq_summary <- rbind(deltarq_summary,data.frame("LABEL"=lbl, tempDataFrame))
  #       else
  #         deltarq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
  #     }
  #     rm(localFileRes)
  #   }
  # }

  if(!plyr::empty(deltarq_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltarq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
