#
#
# library("rjson")
# library(plyr)
#
# json_file <- "~/Downloads/ewas_datahub_metadata-3.txt"
#
# json_data <- fromJSON(paste(readLines(json_file), collapse = ""))
#
# mydata <- data.frame("project id"="")
# mydata <- mydata[mydata$project.id !=NULL]
# for (i in 1:30) {
#   temp <-
#     json_data[which(sapply(json_data, function(x) {
#       length(x) == i
#     }))]
#   if(length(temp)==0)
#     next
#   head_temp <- names(temp[[1]])
#   temp <- matrix(unlist(temp), ncol = i, byrow = TRUE)
#   colnames(temp) <- head_temp
#   mydata <- rbind.fill.matrix(mydata, temp)
# }
#
