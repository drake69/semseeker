#' @importFrom doRNG %dorng%
create_deltaq <- function(envir, resultPopulation){

  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "envir", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
                           ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                           ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
                           ".RNGkind_length", "tail", "RNGstr","update_multiple_bed")
  i <- 2
  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
                           "FIGURE"=envir$keys_figures_default[,1],
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
      deltaq_temp <- fst::read.fst(fileToRead, as.data.table = T)
      colnames(deltaq_temp) <- c("CHR","START","END","VALUE","SAMPLEID")
      deltaq_temp$FIGURE <- key$FIGURE
      deltaq_temp$POPULATION <- key$POPULATION
      if(exists("deltaq"))
        deltaq <- rbind(deltaq, deltaq_temp)
      else
        deltaq <- deltaq_temp
    }
  }

  if (!exists("deltaq")| plyr::empty(deltaq))
  {
     stop("Something wrong with multiple bed files!")
  }

  deltaq$DELTAQ <- as.numeric(dplyr::ntile(x=deltaq[,"VALUE"] , n=4)) * 2

  localKeys <- expand.grid("POPULATION"=unique(resultPopulation$Sample_Group),
                           "FIGURE"=envir$keys_figures_default[,1],
                           "ANOMALY"="DELTAS" ,
                           "EXT"="fst"
  )

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]

    if(exists("deltaq"))
    {
      localFileRes <- deltaq[ deltaq$POPULATION==as.character(key$POPULATION)
                              & deltaq$FIGURE==as.character(key$FIGURE)
                              ,c("CHR","START","END","DELTAQ","SAMPLEID")]
      if(!plyr::empty(localFileRes))
      {
        colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")

        tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste("DELTAQ_",as.character(key$FIGURE),sep="")))
        fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"), as.character(key$FIGURE)),  as.character(key$EXT))
        fst::write.fst( x= localFileRes, fileToWrite)
        message(Sys.time(), " Created DELTAQ multiple annotated file!", fileToWrite)

        tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
        colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
        lbl <- rep(paste("DELTAQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
        data.frame("LABEL"=lbl, tempDataFrame)

        if(exists("deltaq_summary"))
          deltaq_summary <- rbind(deltaq_summary,data.frame("LABEL"=lbl, tempDataFrame))
        else
          deltaq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
      }
    }
  }


  # for(i in 1:nrow(localKeys))
  # {
  #   key <- localKeys[i,]
  #   tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
  #   fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)),  as.character(key$EXT))
  #   if(file.exists(fileToRead))
  #   {
  #     localFileRes <- fst::read.fst(fileToRead, as.data.table = T)
  #
  #     if(!plyr::empty(localFileRes))
  #     {
  #       colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")
  #       deltaq <- unique(deltaq)
  #       localFileRes[,4] <- deltaq[ deltaq$VALUE %in% localFileRes$VALUE, "DELTAQ"]
  #
  #       tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste("DELTAQ_",as.character(key$FIGURE),sep="")))
  #       fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"), as.character(key$FIGURE)),  as.character(key$EXT))
  #       fst::write.fst( x= localFileRes, fileToWrite)
  #       # localFileRes <-fst::read.fst(fileToWrite)
  #
  #       tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
  #       colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
  #       lbl <- rep(paste("DELTAQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
  #       data.frame("LABEL"=lbl, tempDataFrame)
  #
  #       if(exists("deltaq_summary"))
  #         deltaq_summary <- rbind(deltaq_summary,data.frame("LABEL"=lbl, tempDataFrame))
  #       else
  #         deltaq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
  #     }
  #     rm(localFileRes)
  #   }
  # }

  if(!plyr::empty(deltaq_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltaq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
