probe_features_get <- function(area_subarea)
{
  ssEnv <- get_session_info()

  if(is.null(ssEnv$tech) | ssEnv$tech=="")
    stop("ERROR: ", Sys.time(), " Probes get should be called once defined which tech is used.")

  probe_features <- semseeker::pp_tot
  probe_features <- probe_features[ probe_features[,ssEnv$tech], ]
  if(any(area_subarea %in% c("CHR","START")))
    probe_features <- na.omit(probe_features[ c(ssEnv$tech,"PROBE","CHR","START")])
  else
    probe_features <- na.omit(probe_features[ c(ssEnv$tech,"PROBE","CHR","START",area_subarea)])


  probe_features <- probe_features[ ,-c(which(colnames(probe_features) %in% ssEnv$tech )) ]
  probe_features$END <- probe_features$START

  probe_features <- unique(probe_features[,c("PROBE", "CHR","START","END",area_subarea)])
  # colnames(probe_features) <- c("CHR","START","END","AREA")

  return(probe_features)
}
