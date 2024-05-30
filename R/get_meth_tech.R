get_meth_tech <- function(signal_data)
{
  ssEnv <- get_session_info()

  if(nrow(signal_data) == 485512)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems a 450k dataset.")

  if(nrow(signal_data) == 27578)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems a 27k dataset.")

  if(nrow(signal_data) == 866562)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems an EPIC dataset.")

  probe_features <- semseeker::pp_tot[,c("PROBE","CHR","K27","K450","K850")]
  signal_data$PROBE <- rownames(signal_data)
  signal_data_check <- merge(signal_data,probe_features, by="PROBE")

  signal_data_check <- stats::na.omit(subset(signal_data_check, "CHR" !=""))

  tech <- colSums(signal_data_check[,c("K27","K450","K850")])

  tech <-  c("K27","K450","K850")[which(tech==max(tech))]
  if (length(c(tech))>1)
    tech <- tech[1]

  msg = switch(
    tech,
    "K27"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 27k dataset."),
    "K450"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 450k dataset."),
    "K850"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 850k dataset.")
  )

  log_event(msg)
  ssEnv$tech <- tech

  # browser()

  # get abs max
  max_data <- max(abs(signal_data[,-(ncol(signal_data))]), na.rm = TRUE)
  if (max_data>1)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The data is not beta values.")
    ssEnv$beta <- FALSE
  }
  ssEnv$probes_count <- nrow(signal_data)

  update_session_info(ssEnv)

  return(ssEnv)
}
