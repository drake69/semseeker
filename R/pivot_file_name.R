pivot_file_name <- function(marker,figure,area,subarea)
{
  ssEnv <- get_session_info()
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  pivot_subfolder <- dir_check_and_create(reportFolder, marker)
  pivot_file_name <- paste0(marker,"_", figure,"_",area,"_",subarea, sep = "")
  pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_file_name,extension =  ".csv" ,add_gz=TRUE)
  return(pivot_file_name)
}
