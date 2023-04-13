get_meth_tech <- function(methylation_data)
{

  ssEnv <- .pkgglobalenv$ssEnv

  if(nrow(methylation_data) == 485512)
    message("INFO: ", Sys.time(), " seems a 450k dataset.")

  if(nrow(methylation_data) == 27578)
    message("INFO: ", Sys.time(), " seems a 27k dataset.")

  if(nrow(methylation_data) == 866562)
    message("INFO: ", Sys.time(), " seems an EPIC dataset.")

  probes <- semseeker::PROBES_CHR_CHR
  methylation_data$PROBE <- rownames(methylation_data)
  methylation_data_check <- merge(methylation_data,probes, by="PROBE")

  k27 <- nrow(unique(subset(methylation_data_check, !k450 & !k850)))>0
  k450 <- nrow(unique(subset(methylation_data_check, k450 & !k850)))>0
  k850 <- nrow(unique(subset(methylation_data_check, !k450 & k850)))>0

  if(k27)
  {
    ssEnv$tech <- "k27"
    message("INFO: ", Sys.time(), " the dataset is a 27k dataset.")
  }

  if(k450)
  {
    ssEnv$tech <- "k450"
    message("INFO: ", Sys.time(), " the dataset is a 450k dataset.")
  }

  if(k850)
  {
    ssEnv$tech <- "k850"
    message("INFO: ", Sys.time(), " the dataset is a 850k dataset.")
  }

  assign("ssEnv", ssEnv, envir=.pkgglobalenv)

  return(ssEnv)
}
