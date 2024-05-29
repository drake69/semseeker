#' @importFrom doRNG %dorng%
deltap_get <- function(resultPopulation){

  ssEnv <- get_session_info()

  # local_deltap <- paste0("DELTAP",ssEnv$epiquantile, sep="")
  #create multiple file bed
  variables_to_export <- c("localKeys", "resultPopulation", "dir_check_and_create", "ssEnv", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
    ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
    ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
    ".RNGkind_length", "tail", "RNGstr","update_multiple_bed","probe_features_get","log_event")

  Sample_Group=as.data.frame(unique(resultPopulation$Sample_Group))
  colnames(Sample_Group) <- "SAMPLE_GROUP"

  # must use keys_markers_figures_default because the selected marker could exclude deltas which is basic for deltap
  localKeys <- reshape::expand.grid.df(as.data.frame(ssEnv$keys_markers_figures_default),Sample_Group)
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
    tempresult_folderData <- dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste(as.character(key$MARKER),"_",as.character(key$FIGURE),sep="")))
    fileToRead <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$MARKER), as.character(key$FIGURE)), as.character(key$EXT))
    if(file.exists(fileToRead))
    {
      deltas_temp <- fst::read.fst(fileToRead, as.data.table = T)
      colnames(deltas_temp) <- c("CHR","START","END","VALUE","SAMPLEID")
      deltas_temp$FIGURE <- key$FIGURE
      deltas_temp$SAMPLE_GROUP <- key$SAMPLE_GROUP
      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  "deltap_get:", key)
      if(exists("deltap"))
        deltap <- rbind(deltap, deltas_temp)
      else
        deltap <- deltas_temp
    }
    if(ssEnv$showprogress)
      progress_bar(sprintf("Collecting bed files."))
  }

  if (!exists("deltap") | plyr::empty(deltap))
  {
    stop("Something wrong with multiple bed files!")
  }


  # browser()
  deltap$VALUE <- as.numeric(deltap$VALUE)
  num_bins <- as.numeric(ssEnv$bins)
  deltap$DELTAP <- cut(deltap$VALUE, breaks=num_bins, labels=FALSE)
  deltap$DELTAP <- as.numeric(deltap$DELTAP)

  # deltap$DELTAP <- as.numeric(dplyr::ntile(x=deltap[,"VALUE"] , n= as.numeric(ssEnv$epiquantile)))
  deltap_both <- deltap
  deltap_both$FIGURE <- "BOTH"
  deltap_both_sum <- deltap
  deltap_both_sum$DELTAP <- ifelse((deltap_both_sum$FIGURE=="HYPO"), -1 * deltap_both_sum$DELTAP ,deltap_both_sum$DELTAP )
  deltap_both_sum$FIGURE <- "BOTHSUM"
  deltap <- rbind(deltap, deltap_both, deltap_both_sum)

  localKeys <-   reshape::expand.grid.df(as.data.frame(ssEnv$keys_markers_figures_default), data.frame("SAMPLE_GROUP"=unique(resultPopulation$Sample_Group)))
  localKeys$EXT <- "fst"
  localKeys <- subset(localKeys, localKeys$MARKER=="DELTAP")


  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))

  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]

    if(exists("deltap"))
    {
      localFileRes <- deltap[ deltap$SAMPLE_GROUP==as.character(key$SAMPLE_GROUP)
        & deltap$FIGURE==as.character(key$FIGURE)
        ,c("CHR","START","END","DELTAP","SAMPLEID")]
      if(!plyr::empty(localFileRes))
      {
        colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")

        tempresult_folderData <- dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste("DELTAP",key$SUFFIX,"_",as.character(key$FIGURE),sep="")))
        fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAP"),key$SUFFIX, as.character(key$FIGURE)),  as.character(key$EXT))
        fst::write.fst( x= localFileRes, fileToWrite, compress = 100)
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Created DELTAP multiple annotated file!", fileToWrite)

        tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",fun.aggregate = sum, drop = TRUE)
        colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
        lbl <- rep(paste("DELTAP", key$SUFFIX,"_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
        data.frame("LABEL"=lbl, tempDataFrame)

        if(exists("deltap_summary"))
          deltap_summary <- rbind(deltap_summary,data.frame("LABEL"=lbl, tempDataFrame))
        else
          deltap_summary <- data.frame("LABEL"=lbl, tempDataFrame)
      }
    }
    if(ssEnv$showprogress)
      progress_bar(sprintf("Binding bed files."))

  }

  if(!plyr::empty(deltap_summary))
  {
    tempDataFrame <- reshape2::dcast(data = deltap_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)
    # remove from result columns existing in tempDataFrame
    resultPopulation <- resultPopulation[,!(colnames(resultPopulation) %in% colnames(tempDataFrame))]
    resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")
  }
  return(resultPopulation)

}
