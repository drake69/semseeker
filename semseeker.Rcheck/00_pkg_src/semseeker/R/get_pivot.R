get_pivot <- function(marker, figure, group, subgroup){

  ssEnv <- get_session_info()

  pivot_folder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  pivot_file_name <-  paste(marker,"_",figure,"_",group,"_",subgroup, sep="")
  fileName <- paste0(pivot_folder,"/",pivot_file_name,".csv" , sep="")

  pivot_data <- utils::read.csv2(fileName, sep=";")

  return(pivot_data)
}
