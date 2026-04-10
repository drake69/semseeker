probe_features_get <- function(area_subarea)
{
  ssEnv <- get_session_info()

  if(is.null(ssEnv$tech) | ssEnv$tech=="")
  {

    # try to get from signal_data saved
    pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE","WHOLE")
    signal_data_pl <- polars::pl$read_parquet(pivot_file_name)
    signal_data_r  <- as.data.frame(signal_data_pl)
    # probe names are stored in the AREA column; restore them as rownames
    if ("AREA" %in% colnames(signal_data_r)) {
      rownames(signal_data_r) <- signal_data_r$AREA
      signal_data_r$AREA <- NULL
    }
    ssEnv <- get_meth_tech(signal_data_r)
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Probes get should be called once defined which tech is used.")
    if(is.null(ssEnv$tech) | ssEnv$tech=="")
    {
      log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Probes get should be called once defined which tech is used.")
      stop()
    }
  }


  if(!grepl("_",area_subarea))
    area_subarea <- paste(area_subarea,"_","WHOLE",sep="")

  probe_features <- SEMseeker::pp_tot
  probe_features <- probe_features[ probe_features[,ssEnv$tech], ]
  probe_features$END <- probe_features$START

  if( (grepl("CHR", area_subarea)  || grepl("PROBE",area_subarea)) && !grepl("CHR_CYTOBAND",area_subarea))
    # probe_features <- unique(probe_features[ c(ssEnv$tech,"PROBE","CHR","START","END")])
    probe_features <- probe_features[ c(ssEnv$tech,"PROBE","CHR","START","END")] %>% dplyr::distinct()
  else
    # probe_features <- unique(probe_features[ c(ssEnv$tech,"PROBE","CHR","START","END",area_subarea)])
    probe_features <- probe_features[ c(ssEnv$tech,"PROBE","CHR","START","END",area_subarea)] %>% dplyr::distinct()

  # needed to avoid special elaboration in other function on chromosome or probe
  if( grepl("CHR", area_subarea) && !grepl("CHR_CYTOBAND",area_subarea))
    probe_features$CHR_WHOLE <- paste("chr",probe_features$CHR, sep="")

  if(grepl("PROBE",area_subarea))
    probe_features$PROBE_WHOLE <- probe_features$PROBE

  probe_features <- probe_features[ ,-c(which(colnames(probe_features) %in% ssEnv$tech )) ]

  # remove sex chromosomes
  if(ssEnv$sex_chromosome_remove==TRUE)
    probe_features <- probe_features[ !(probe_features$CHR %in% c("X","Y")), ]

  return(probe_features)
}
