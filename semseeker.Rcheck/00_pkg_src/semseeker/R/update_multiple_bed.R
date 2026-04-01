# update_multiple_bed <- function (path_file_name, data_to_add)
# {
#   if(file.exists(path_file_name))
#   {
#     data <- fst::read.fst(path_file_name, as.data.table = T)
#     data <- rbind(data, data_to_add)
#   } else
#     data <- data_to_add
#
#   fst::write.fst(x = data,path = path_file_name)
#
# }
# #
