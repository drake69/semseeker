#' @importFrom doRNG %dorng%
deltaq_get <- function(resultPopulation){

  ssEnv <- semseeker:::get_session_info()

  # local_deltaq <- paste0("DELTAQ",ssEnv$epiquantile, sep="")
  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "ssEnv", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
                           ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                           ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
                           ".RNGkind_length", "tail", "RNGstr","update_multiple_bed","semseeker:::probe_features_get")

  Sample_Group=as.data.frame(unique(resultPopulation$Sample_Group))
  colnames(Sample_Group) <- "SAMPLE_GROUP"

  # must use keys_markers_figure_default because the selected marker could exclude deltas which is basic for deltaq
  localKeys <- reshape::expand.grid.df(ssEnv$keys_markers_figure_default,Sample_Group)
  localKeys <- subset(localKeys, localKeys$MARKER=="DELTAS")
  localKeys <- subset(localKeys, localKeys$FIGURE!="BOTHSUM")
  localKeys <- subset(localKeys, localKeys$FIGURE!="BOTH")
  localKeys$EXT <- "fst"

  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste(as.character(key$MARKER),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- semseeker:::file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$MARKER), as.character(key$FIGURE)), as.character(key$EXT))
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
    if(ssEnv$showprogress)
      progress_bar(sprintf("Collecting bed files."))
  }

  if (!exists("deltaq") | plyr::empty(deltaq))
  {
     stop("Something wrong with multiple bed files!")
  }

  deltaq$DELTAQ <- as.numeric(dplyr::ntile(x=deltaq[,"VALUE"] , n= as.numeric(ssEnv$epiquantile)))
  deltaq_both <- deltaq
  deltaq_both$FIGURE <- "BOTH"
  deltaq_both_sum <- deltaq
  deltaq_both_sum$DELTAQ <- ifelse((deltaq_both_sum$FIGURE=="HYPO"), -1 * deltaq_both_sum$DELTAQ ,deltaq_both_sum$DELTAQ )
  deltaq_both_sum$FIGURE <- "BOTHSUM"
  deltaq <- rbind(deltaq, deltaq_both, deltaq_both_sum)

  localKeys <-   reshape::expand.grid.df(ssEnv$keys_markers_figure_default, data.frame("SAMPLE_GROUP"=unique(resultPopulation$Sample_Group)))
  localKeys$EXT <- "fst"
  localKeys <- subset(localKeys, localKeys$MARKER=="DELTAQ")

  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))

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

        tempresult_folderData <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste("DELTAQ",key$SUFFIX,"_",as.character(key$FIGURE),sep="")))
        fileToWrite <- semseeker:::file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"),key$SUFFIX, as.character(key$FIGURE)),  as.character(key$EXT))
        fst::write.fst( x= localFileRes, fileToWrite, compress = 100)
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Created DELTAQ multiple annotated file!", fileToWrite)

        tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
        colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
        lbl <- rep(paste("DELTAQ", key$SUFFIX,"_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
        data.frame("LABEL"=lbl, tempDataFrame)

        if(exists("deltaq_summary"))
          deltaq_summary <- rbind(deltaq_summary,data.frame("LABEL"=lbl, tempDataFrame))
        else
          deltaq_summary <- data.frame("LABEL"=lbl, tempDataFrame)
      }
    }
    if(ssEnv$showprogress)
      progress_bar(sprintf("Binding bed files."))

  }

  if(!plyr::empty(deltaq_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltaq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    # remove from result columns existing in tempDataFrame
    resultPopulation <- resultPopulation[,!(colnames(resultPopulation) %in% colnames(tempDataFrame))]
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
