#' @importFrom doRNG %dorng%
create_multiple_bed <- function(sample_sheet){

  ssEnv <- get_session_info()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Started multiple file creation!")
  #create multiple file bed
  i <- 0
  Sample_Group <- as.data.frame(unique(sample_sheet$Sample_Group))
  colnames(Sample_Group) <- "SAMPLE_GROUP"
  # multiple bed i created only for probes must be taken into account all markers
  # 
  localKeys <- as.data.frame(reshape::expand.grid.df(as.data.frame(ssEnv$keys_markers_figures_default),Sample_Group))
  # localKeys <- subset(localKeys, MARKER!="SIGNAL")

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(localKeys))
  else
    progress_bar <- ""

  to_export <- c("localKeys", "dir_check_and_create", "ssEnv", "file_path_build", "sample_sheet","progress_bar",
    "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","log_event")

  foreach::foreach(i = 1:nrow(localKeys), .export = to_export) %dorng%
  # for(i in 1:nrow(localKeys))
  {
    # 
    # i <- 1
    key <- localKeys[i,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    tempresult_folderData <-dir_check_and_create(ssEnv$result_folderData,c(as.character(key$SAMPLE_GROUP) ,paste(as.character(marker),"_",as.character(figure),sep="")))
    temp_file <- tempdir()
    temp_file <- paste(temp_file, stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),sep="")
    fileToWrite <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(marker), as.character(figure)), "fst")
    fileToWriteBed <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(marker), key$SUFFIX, as.character(figure)), key$EXT, add_gz=TRUE)
    if(!file.exists(fileToWrite))
    {
      j <- 0
      for ( j in 1:nrow(sample_sheet))
      {
        # j <- 1
        sample <- sample_sheet[j,]
        fileToRead <- file_path_build(baseFolder =  tempresult_folderData,
          detailsFilename =  c(sample$Sample_ID, as.character(marker),key$SUFFIX, as.character(figure)),
          extension =  key$EXT,
          add_gz = TRUE)
        if(file.exists(fileToRead))
        {
          localtemp <- utils::read.table(gzfile(fileToRead), sep="\t", header = FALSE)
          # localtemp <- utils::read.csv2(fileToRead, sep="\t", header = FALSE)
          localtemp$Sample_ID <- sample$Sample_ID
          if(!plyr::empty(localtemp))
          {
            utils::write.table(localtemp, temp_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
          }
          rm(localtemp)
        }
      }

      # 
      if(file.exists(temp_file))
      {
        fst::write.fst(x = utils::read.table(temp_file, sep = "\t"),path = fileToWrite, compress = 100)
        # read temp file
        localtemp <- utils::read.table(temp_file, sep="\t", header = FALSE)
        # write gz file
        utils::write.table(localtemp,gzfile(fileToWriteBed), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
        file.remove(temp_file)
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Created multiple annotated file!", fileToWriteBed)
      }
    }
    if(ssEnv$showprogress)
      progress_bar(sprintf("Creating multiple file."))
  }
  # future::plan(ssEnv$parallel_strategy)
}
