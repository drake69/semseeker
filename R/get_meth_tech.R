get_meth_tech <- function(signal_data)
{
  ssEnv <- get_session_info()

  if(nrow(signal_data) == 485512)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems a 450k dataset.")

  if(nrow(signal_data) == 27578)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems a 27k dataset.")

  if(nrow(signal_data) == 866562)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems an EPIC dataset.")

  if(nrow(signal_data) > 866562)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " seems a WGBS dataset.")

  # if doesn't exist column PROBE
  if(!"PROBE" %in% colnames(signal_data))
    signal_data_probe <- rownames(signal_data)
  else
    signal_data_probe <- signal_data$PROBE

  probe_features <- semseeker::pp_tot[,c("PROBE","CHR","K27","K450","K850")]
  # get only probe_features that are in the signal_data_probe
  probe_features <- probe_features[probe_features$PROBE %in% signal_data_probe,]
  # signal_data_check <- merge(signal_data,probe_features, by="PROBE")
  # if(exists("signal_data"))
  #   rm(signal_data)
  # rm(probe_features)

  tech <- ""
  msg <- ""
  if (nrow(probe_features)!=0)
  {
    probe_features <- subset(probe_features, "CHR" !="")
    tech <- colSums(probe_features[,c("K27","K450","K850")])
    tech <-  c("K27","K450","K850")[which(tech==max(tech))]
    if (length(c(tech))>=1)
    {
      tech <- tech[1]
      msg <- switch(
        tech,
        "K27"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 27k dataset."),
        "K450"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 450k dataset."),
        "K850"= paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a 850k dataset.")
      )
    }
  }

  if (all(grepl("_", probe_features$PROBE[1:1000])))
  {
    tech <- "WGBS"
    msg <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is a WGBS dataset.")
  }

  if(tech == "")
  {
    tech <- "UNKNOWN"
    msg <- paste("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), "the dataset is an unknown dataset.")
    log_event(msg)
    stop(msg)
  }

  log_event(msg)
  ssEnv$tech <- tech
  # remove column PROBE
  signal_data_check <- signal_data
  if("PROBE" %in% colnames(signal_data_check))
    signal_data_check <- signal_data_check[,-which(colnames(signal_data_check)=="PROBE")]

  log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " defining if meth is a beta or mvalue.")
  exploratory_df <- signal_data_check[1:min(10000,nrow(signal_data_check)),!(colnames(signal_data_check) %in% c("PROBE","CHR","K27","K450","K850"))]
  # get abs max
  max_data <- max(exploratory_df, na.rm=TRUE)
  min_data <- min(exploratory_df, na.rm=TRUE)
  max_data <- max(abs(c(max_data, min_data)))
  ssEnv$beta <- TRUE
  if (max_data>1)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The data are mValues.")
    ssEnv$beta <- FALSE
  }
  else
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " The data are beta values.")
  ssEnv$probes_count <- nrow(signal_data_check)

  update_session_info(ssEnv)

  return(ssEnv)
}
