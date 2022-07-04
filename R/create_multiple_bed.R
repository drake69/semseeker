#' @importFrom doRNG %dorng%
create_multiple_bed <- function(envir, sample_sheet, resultPopulation){

  #create multiple file bed
  variables_to_export <- c("localKeys", "sample_sheet", "dir_check_and_create", "envir", "file_path_build","%dorng%","getdorng","iter", "RNGseed", "checkRNGversion", "getRNG", "%||%",
                           ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                           ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider",
                           ".RNGkind_length", "tail", "RNGstr")
  i <- 0



  localKeys <- expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),"FIGURE"=envir$keys_figures_default[,1] ,"ANOMALY"= envir$keys_anomalies_default[,1] ,"EXT"="bed")
  localKeys <- rbind(localKeys, expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
                                            "FIGURE"=envir$keys_figures_default[,1],
                                            "ANOMALY"="DELTAS" ,"EXT"="bedgraph"))

  deltaq_summary <- foreach::foreach(i = 1:nrow(localKeys), .export = variables_to_export, .combine = rbind ) %dorng%
  # for(i in 1:nrow(localKeys))
    {
      key <- localKeys[i,]
      tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
      variables_to_export_nested <- c("sample_sheet", "file_path_build","key","tempresult_folderData","%dorng%","getdorng")
      j <- 0
      localFileRes <- foreach::foreach(j = 1:nrow(sample_sheet), .export = variables_to_export_nested, .combine = "rbind" ) %dorng%
        {
          sample <- sample_sheet[j,]
          fileToRead <- file_path_build(tempresult_folderData, c(sample$Sample_ID, as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
          if(file.exists(fileToRead))
          {
            localtemp <- utils::read.csv2(fileToRead, sep="\t", header = FALSE)
            localtemp$Sample_ID <- sample$Sample_ID
            localtemp
          }
        }

      if(exists("localFileRes"))
        if(!plyr::empty(localFileRes))
        {
          fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
          utils::write.table(localFileRes, fileToWrite, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
          #create quantilized deltas
          if(key$ANOMALY=="DELTAS" & sum("DELTAQ" %in% envir$keys_anomalies)>0)
          {
            #value is in the 4th position
            # give to each quantile an even weight
            localFileRes[,4] <- as.numeric(dplyr::ntile(x=as.numeric(localFileRes[,4]) , n=4)) * 2
            tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste("DELTAQ_",as.character(key$FIGURE),sep="")))
            fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character("DELTAQ"), as.character(key$FIGURE)), key$EXT)
            utils::write.table(localFileRes, fileToWrite, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
            colnames(localFileRes) <- c("CHR","START","END","VALUE","SAMPLEID")
            tempDataFrame <- reshape2::dcast(data = localFileRes, formula = SAMPLEID  ~ ., value.var = "VALUE",
                                             fun.aggregate = sum, drop = TRUE)
            colnames(tempDataFrame) <- c("SAMPLEID","VALUE")
            lbl <- rep(paste("DELTAQ_",as.character(key$FIGURE),sep=""), nrow(tempDataFrame))
            data.frame("LABEL"=lbl, tempDataFrame)
          }
        }
    }

  #update sample sheet
  # study_summary <-   utils::read.csv2(file_path_build( envir$result_folderData, "sample_sheet_result","csv"))

  tempDataFrame <- reshape2::dcast(data = deltaq_summary, formula = SAMPLEID  ~ LABEL, value.var = "VALUE", fun.aggregate = sum, drop = TRUE)

  # study_summary <- study_summary[, !(colnames(study_summary) %in% colnames(tempDataFrame))]
  resultPopulation <- merge(resultPopulation, tempDataFrame, by.x="Sample_ID", by.y="SAMPLEID")

  gc()
  return(resultPopulation)

}
