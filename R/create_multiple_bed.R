#' @importFrom doRNG %dorng%
create_multiple_bed <- function(sample_sheet){

  ssEnv <- get_session_info()

  message("INFO: ", Sys.time(), " Started multiple annotated file creation!")

  #create multiple file bed
  i <- 0
  localKeys <- expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
                           "FIGURE"=ssEnv$keys_figures_default[,1] ,
                           "ANOMALY"= c("MUTATIONS","LESIONS","DELTAQ","DELTARQ") ,
                           "EXT"="bed")
  localKeys <- rbind(localKeys, expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
                                            "FIGURE"=ssEnv$keys_figures_default[,1],
                                            "ANOMALY"=c("DELTAS","DELTAR") ,
                                            "EXT"="bedgraph")
                    )
  # localKeys <- rbind(localKeys, expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
  #   "FIGURE"=ssEnv$keys_figures_default[,1],
  #   "ANOMALY"=c("DELTAS","DELTAR","BETA") ,
  #   "EXT"="bedgraph")
  # )

  to_export <- c("localKeys", "dir_check_and_create", "ssEnv", "file_path_build", "sample_sheet")
  # future::plan( future::sequential)
  # foreach::foreach(i = 1:nrow(localKeys), .export = to_export) %dorng%
  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    anomaly <- as.character(key$ANOMALY)
    figure <- as.character(key$FIGURE)
    if (anomaly=="BETA")
    {
      figure <- "MEAN"
    }
    tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(anomaly),"_",as.character(figure),sep="")))
    temp_file <- tempdir()
    temp_file <- paste(temp_file, stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),sep="")
    fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(anomaly), as.character(figure)), "fst")
    if(!file.exists(fileToWrite))
    {
      j <- 0
      for ( j in 1:nrow(sample_sheet))
      {
        sample <- sample_sheet[j,]
        fileToRead <- file_path_build(tempresult_folderData, c(sample$Sample_ID, as.character(anomaly), as.character(figure)), key$EXT)
        if(file.exists(fileToRead))
        {
          localtemp <- utils::read.csv2(fileToRead, sep="\t", header = FALSE)
          localtemp$Sample_ID <- sample$Sample_ID
          if(!plyr::empty(localtemp))
          {
            utils::write.table(localtemp, temp_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
          }
          rm(localtemp)
        }
      }

      if(file.exists(temp_file))
      {
        # if(ncol(temp_file)==4)
        #   colnames(temp_file)<- c("CHR","START","END","SAMPLEID")
        # else
        #   colnames(temp_file)<- c("CHR","START","END","VALUE","SAMPLEID")

        fst::write.fst(x = utils::read.table(temp_file, sep = "\t"),path = fileToWrite)
        file.remove(temp_file)
        message("INFO: ", Sys.time(), " Created multiple annotated file!", fileToWrite)
      }
    }
  }
  # future::plan(ssEnv$parallel_strategy)
}
