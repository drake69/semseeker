#' @importFrom doRNG %dorng%
deltaq_get <- function(resultPopulation){

  ssEnv <- get_session_info()

  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "ssEnv", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
                           ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                           ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
                           ".RNGkind_length", "tail", "RNGstr","update_multiple_bed","probe_features_get")

  Sample_Group=as.data.frame(unique(resultPopulation$Sample_Group))
  colnames(Sample_Group) <- "SAMPLE_GROUP"
  localKeys <- reshape::expand.grid.df(ssEnv$keys_markers_figures,Sample_Group)
  localKeys <- subset(localKeys, localKeys$MARKER=="DELTAS")
  localKeys$EXT <- "fst"

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste(as.character(key$MARKER),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$MARKER), as.character(key$FIGURE)), as.character(key$EXT))
    if(file.exists(fileToRead))
    {
      deltaq_temp <- fst::read.fst(fileToRead, as.data.table = T)
      colnames(deltaq_temp) <- c("CHR","START","END","VALUE","SAMPLEID")
      deltaq_temp$FIGURE <- key$FIGURE
      deltaq_temp$SAMPLE_GROUP <- key$SAMPLE_GROUP
      if(exists("deltaq"))
        deltaq <- rbind(deltaq, deltaq_temp)
      else
        deltaq <- deltaq_temp
    }
  }

  if (!exists("deltaq") | plyr::empty(deltaq))
  {
     stop("Something wrong with multiple bed files!")
  }

  deltaq$DELTAQ <- as.numeric(dplyr::ntile(x=deltaq[,"VALUE"] , n=4))
  localKeys <-   reshape::expand.grid.df(ssEnv$keys_markers_figures, data.frame("SAMPLE_GROUP"=unique(resultPopulation$Sample_Group)))
  localKeys$EXT <- "fst"
  localKeys <- subset(localKeys, localKeys$MARKER=="DELTAQ")

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]

    if(exists("deltaq"))
    {
      localFileRes <- deltaq[ deltaq$SAMPLE_GROUP==as.character(key$SAMPLE_GROUP)
                              & deltaq$FIGURE==as.character(key$FIGURE)
                              ,c("CHR","START","END","DELTAQ","SAMPLEID")]
      if(!plyr::empty(localFileRes))
      {
        colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")

        tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste("DELTAQ_",as.character(key$FIGURE),sep="")))
        fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"), as.character(key$FIGURE)),  as.character(key$EXT))
        fst::write.fst( x= localFileRes, fileToWrite)
        message("INFO: ", Sys.time(), " Created DELTAQ multiple annotated file!", fileToWrite)

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

  if(!plyr::empty(deltaq_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltaq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
