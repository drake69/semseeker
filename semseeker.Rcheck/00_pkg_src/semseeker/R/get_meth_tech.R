get_meth_tech <- function(methylation_data)
{
  ssEnv <- get_session_info()

  if(nrow(methylation_data) == 485512)
    message("INFO: ", Sys.time(), " seems a 450k dataset.")

  if(nrow(methylation_data) == 27578)
    message("INFO: ", Sys.time(), " seems a 27k dataset.")

  if(nrow(methylation_data) == 866562)
    message("INFO: ", Sys.time(), " seems an EPIC dataset.")

  probe_features <- semseeker::pp_tot[,c("PROBE","CHR","K27","K450","K850")]
  methylation_data$PROBE <- rownames(methylation_data)
  methylation_data_check <- merge(methylation_data,probe_features, by="PROBE")

  methylation_data_check <- stats::na.omit(subset(methylation_data_check, "CHR" !=""))

  tech <- colSums(methylation_data_check[,c("K27","K450","K850")])

  tech <-  c("K27","K450","K850")[which(tech==max(tech))]
  if (length(c(tech))>1)
    tech <- tech[1]

  msg = switch(
    tech,
    "K27"= paste("INFO:", Sys.time(), "the dataset is a 27k dataset."),
    "K450"= paste("INFO:", Sys.time(), "the dataset is a 450k dataset."),
    "K850"= paste("INFO:", Sys.time(), "the dataset is a 850k dataset.")
  )

  message(msg)
  ssEnv$tech <- tech
  update_session_info(ssEnv)

  return(ssEnv)
}
