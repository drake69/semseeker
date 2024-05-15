analyze_single_sample_both <- function( sample_detail, marker) {

  ssEnv <- get_session_info()
  # log_event(sample_detail$Sample_ID, " ", "SingleSampleBoth Sample analysis warmingUP ", format(Sys.time(), "%a %b %d %X %Y"))
  result <- data.frame()

  # nothing to do
  if ( !any(ssEnv$keys_markers_figures$COMBINED==paste(marker,"_","BOTH",sep="")) &
      !any(ssEnv$keys_markers_figures$COMBINED==paste(marker,"_","BOTHSUM",sep="")))
    return(result)

  start_time_single_sample <- Sys.time()

  data_to_save <- data.frame("CHR"="", "START"="", "END"="", "VALUE"="")
  data_to_save <- data_to_save[-1,]
  marker <- as.character(marker)
  figures <- c("HYPER","HYPO")

  for( i in 1:length(figures))
  {
    figure <- figures[i]
    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0(marker,"_", figure, sep = "")))
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,marker,figure),"bed", add_gz=TRUE)
    if(file.exists(fileName))
    {
      data_to_saveTemp <- utils::read.table(fileName, sep="\t", col.names =c("CHR", "START", "END","VALUE") )
      if (figure == "HYPER")
        data_to_saveTemp$VALUE <- 1
      else
        data_to_saveTemp$VALUE <- -1
      data_to_save <- plyr::rbind.fill(data_to_save, data_to_saveTemp )
    }
  }

  if(nrow(data_to_save) == 0)
  {
    end_time_single_sample <- Sys.time()
    time_taken <- difftime(end_time_single_sample,start_time_single_sample, units = "mins")
    log_event("DEBUG: Completed figure BOTH for sample in ", as.numeric(time_taken), "minutes.")
    return()
  }

  if (any(ssEnv$keys_markers_figures$COMBINED==paste(marker,"_","BOTH",sep="")))
  {
    figure <- "BOTH"
    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0(marker,"_", figure, sep = "")))
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,marker,figure),"bed", add_gz=TRUE)
    data_to_save_both <- data_to_save
    data_to_save_both$VALUE <- 1

    dump_sample_as_bed_file(
      data_to_dump = data_to_save,
      fileName = fileName
    )
    result <- data.frame_add.column(result, paste(marker,"_","BOTH",sep=""),sum(data_to_save_both$VALUE))
  }

  if (any(ssEnv$keys_markers_figures$COMBINED==paste(marker,"_","BOTHSUM",sep="")))
  {
    figure <- "BOTHSUM"
    folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0(marker,"_", figure, sep = "")))
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,marker,figure),"bed", add_gz=TRUE)
    dump_sample_as_bed_file(
      data_to_dump = data_to_save,
      fileName = fileName
    )
    result <- data.frame_add.column(result, paste(marker,"_","BOTHSUM",sep=""),sum(as.numeric(data_to_save$VALUE)))

  }

  end_time_single_sample <- Sys.time()
  time_taken <- difftime(end_time_single_sample,start_time_single_sample, units = "mins")
  log_event("DEBUG: Completed figure BOTH for sample in ", as.numeric(time_taken), "minutes.")

  return(result)
}
