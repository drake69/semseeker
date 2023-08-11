beta_save <- function(methylation_data, sample_sheet, batch_id )
{
  ssEnv <- get_session_info()


  # save beta as rds and as pivot
  sample_sheet <- subset(sample_sheet, sample_sheet$Sample_Group != "Reference")
  sample_info <- sample_sheet[,c("Sample_ID","Sample_Group")]
  colnames(sample_info) <- c("SAMPLEID","SAMPLE_GROUP")
  methylation_data <- data.frame("SAMPLEID"=rownames(methylation_data), methylation_data[, sample_info$SAMPLEID])
  sample_info <- as.data.frame(t(sample_info))
  colnames(sample_info) <- sample_info[1,]
  sample_info <- cbind(data.frame("SAMPLEID"="SAMPLE_GROUP"), sample_info[-1,])
  methylation_data <- rbind(sample_info, methylation_data )
  pivot_subfolder <- dir_check_and_create(ssEnv$result_folderData,c("Pivots","BETA"))
  fileName <- file_path_build(pivot_subfolder,c("BETA","MEAN","PROBE","PROBE"),"csv")
  utils::write.table(methylation_data, fileName, row.names = T, col.names = T, sep=";")
  saveRDS(methylation_data,file_path_build(ssEnv$result_folderData, c(batch_id,"_methylation_data"),"rds"))

}
