probes_get <- function(probes_prefix, grp)
{

  ssEnv <- get_session_info()

  if(ssEnv$tech=="")
    stop("ERROR: ", Sys.time(), " Probes get should be called once defined which tech is used.")

  probes_filter <- semseeker::PROBES
  probes_filter <- probes_filter[probes_filter[,ssEnv$tech],"PROBE"]

  if(probes_prefix=="PROBES_CHR_")
  {
    probes <- semseeker::PROBES_CHR_CHR
    probes <- probes[probes[,ssEnv$tech],]
    # message("DEBUG: loaded probes: PROBES_CHR_CHR")
  }
  else if(probes_prefix=="PROBES")
  {
    probes <- semseeker::PROBES
    probes <- probes[probes[,ssEnv$tech],]
    # message("DEBUG: loaded probes: PROBES_CHR_CHR")
  }
  else
  {
    probes_name <- paste0(probes_prefix, grp,sep="")
    probes <- get(probes_name, envir = asNamespace("semseeker"))
    # message("DEBUG: loaded probes:", probes_name)
  }

  probes <- unique(probes[probes$PROBE %in% probes_filter, !(colnames(probes) %in% c("ACCESSION","POSITION"))])


  # is ok to have more instance of the same probe because a probe can fall over more than one gene
  # and also a gene can have different isoform
  # if(nrow(probes) != length(unique(probes$PROBE)))
  #   stop("ERROR: ", Sys.Date(), " Probes error!")

  return(probes)
}
