probe_features_get <- function(area_subarea)
{
  ssEnv <- get_session_info()

  if(is.null(ssEnv$tech) | ssEnv$tech=="")
    {
      log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Probes get should be called once defined which tech is used.")
      stop()
    }
  probe_features <- semseeker::pp_tot
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
