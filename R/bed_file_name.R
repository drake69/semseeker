bed_file_name <- function(sample_id,sample_group,marker,figure)
{
  ssEnv <- get_session_info()

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_group),paste0(marker,"_",figure)))
  bed_ext <- unique(ssEnv$keys_markers_figures[ ssEnv$keys_markers_figures$MARKER==marker & ssEnv$keys_markers_figures$FIGURE==figure, "EXT"])
  if(length(bed_ext)==0){
    return()
  }
  file_name <- file_path_build(folder_to_save,c(sample_id,marker,figure),bed_ext, add_gz=TRUE)
  return(file_name)
}
