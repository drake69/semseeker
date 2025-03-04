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

  signal_data$PROBE <- rownames(signal_data)

  probe_features <- semseeker::pp_tot[,c("PROBE","CHR","K27","K450","K850")]
  signal_data_check <- merge(signal_data,probe_features, by="PROBE")
  rm(signal_data)
  rm(probe_features)

  tech <- ""
  msg <- ""
  if (nrow(signal_data_check)!=0)
  {
    signal_data_check <- subset(signal_data_check, "CHR" !="")
    tech <- colSums(signal_data_check[,c("K27","K450","K850")])
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

  if (all(grepl("_", signal_data_check$PROBE[1:1000])))
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

  log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " defining if meth is a beta or mvalue.")
  exploratory_df <- signal_data_check[1:min(10000,nrow(signal_data_check)),!(colnames(signal_data_check) %in% c("PROBE","CHR","K27","K450","K850"))]
  # get abs max
  max_data <- max(exploratory_df)
  min_data <- min(exploratory_df)
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
