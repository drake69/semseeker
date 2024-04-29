signal_save <- function(signal_data, sample_sheet, batch_id )
{
  ssEnv <- semseeker:::get_session_info()

  # save signal as rds and as pivot
  sample_sheet <- subset(sample_sheet, sample_sheet$Sample_Group != "Reference")
  sample_info <- sample_sheet[,c("Sample_ID","Sample_Group")]
  colnames(sample_info) <- c("SAMPLEID","SAMPLE_GROUP")
  signal_data <- data.frame("SAMPLEID"=rownames(signal_data), signal_data[, sample_info$SAMPLEID])
  sample_info <- as.data.frame(t(sample_info))
  colnames(sample_info) <- sample_info[1,]
  sample_info <- cbind(data.frame("SAMPLEID"="SAMPLE_GROUP"), sample_info[-1,])
  signal_data <- signal_data[,colnames(sample_info)]
  signal_data <- plyr::rbind.fill(sample_info, signal_data )
  pivot_subfolder <- semseeker:::dir_check_and_create(ssEnv$result_folderData,c("Pivots","SIGNAL"))
  fileName <- semseeker:::file_path_build(pivot_subfolder,c("SIGNAL","MEAN","PROBE","PROBE"),"csv", add_gz = TRUE)
  utils::write.table(signal_data, gzfile(fileName), row.names = T, col.names = T, sep=";")
  saveRDS(signal_data,file_path_build(ssEnv$result_folderData, c(batch_id,"_signal_data"),"rds"))

}
