get_meth_tech <- function(methylation_data)
{
  ssEnv <- .pkgglobalenv$ssEnv

  if(nrow(methylation_data) == 485512)
    message("INFO: ", Sys.time(), " seems a 450k dataset.")

  if(nrow(methylation_data) == 27578)
    message("INFO: ", Sys.time(), " seems a 27k dataset.")

  if(nrow(methylation_data) == 866562)
    message("INFO: ", Sys.time(), " seems an EPIC dataset.")

  probes <- semseeker::PROBES
  methylation_data$PROBE <- rownames(methylation_data)
  methylation_data_check <- merge(methylation_data,probes, by="PROBE")

  methylation_data_check <- stats::na.omit(subset(methylation_data_check, "CHR" !=""))

  tech <- colSums(methylation_data_check[,c("k27","k450","k850")])

  tech <-  c("k27","k450","k850")[which(tech==max(tech))]
  if (length(c(tech))>1)
    tech <- tech[1]

  msg = switch(
    tech,
    "k27"= paste("INFO: ", Sys.time(), " the dataset is a 27k dataset."),
    "k450"= paste("INFO: ", Sys.time(), " the dataset is a 450k dataset."),
    "k850"= paste("INFO: ", Sys.time(), " the dataset is a 850k dataset.")
  )

  message(msg)
  ssEnv$tech <- tech
  update_session_info(ssEnv)

  return(ssEnv)
}
