#' @importFrom foreach %dopar%
create_multiple_bed <- function(envir, sample_sheet){

  #create multiple file bed
  variables_to_export <- c("localKeys", "sample_sheet", "dir_check_and_create", "envir", "file_path_build","%dopar%","getDoPar")
  i <- 0
  localKeys <- expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),"FIGURE"=envir$keys_figures[,1] ,"ANOMALY"= envir$keys_anomalies[,1] ,"EXT"="bed")
  localKeys <- rbind(localKeys, expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),"FIGURE"=c("HYPO","HYPER","BOTH"),"ANOMALY"="DELTAS" ,"EXT"="bedgraph"))

  foreach::foreach(i = 1:nrow(localKeys), .export = variables_to_export ) %dopar%
    {
      key <- localKeys[i,]
      tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
      variables_to_export_nested <- c("sample_sheet", "file_path_build","key","tempresult_folderData","%dopar%","getDoPar")
      j <- 0
      localFileRes <- foreach::foreach(j = 1:nrow(sample_sheet), .export = variables_to_export_nested, .combine = "rbind" ) %dopar%
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
        if(!plyr::empty(localFileRes)>0)
        {
          fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
          utils::write.table(localFileRes, fileToWrite, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
    }
  gc()
}
