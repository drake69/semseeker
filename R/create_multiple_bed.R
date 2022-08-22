#' @importFrom doRNG %dorng%
create_multiple_bed <- function(envir, sample_sheet){

  message(Sys.time(), " Started multiple annotated file creation!")

  #create multiple file bed
  i <- 0
  localKeys <- expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
                           "FIGURE"=envir$keys_figures_default[,1] ,
                           "ANOMALY"= c("MUTATIONS","LESIONS") ,"EXT"="bed")
  localKeys <- rbind(localKeys, expand.grid("POPULATION"=unique(sample_sheet$Sample_Group),
                                            "FIGURE"=envir$keys_figures_default[,1],
                                            "ANOMALY"="DELTAS" ,"EXT"="bedgraph"))

  # to_export <- c("localKeys", "dir_check_and_create", "envir", "file_path_build", "sample_sheet")
  future::plan( future::sequential)
  # foreach::foreach(i = 1:nrow(localKeys), .export = to_export) %dorng%
  for(i in 1:nrow(localKeys))
  {
    key <- localKeys[i,]
    tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(key$POPULATION) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
    temp_file <- paste("/tmp/", stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),sep="")
    fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), "fst")
    j <- 0
    for ( j in 1:nrow(sample_sheet))
    {
      sample <- sample_sheet[j,]
      fileToRead <- file_path_build(tempresult_folderData, c(sample$Sample_ID, as.character(key$ANOMALY), as.character(key$FIGURE)), key$EXT)
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
        message(Sys.time(), " Created multiple annotated file!", fileToWrite)
    }
    gc()
  }
  future::plan( envir$parallel_strategy)
}
