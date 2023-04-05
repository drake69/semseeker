probes_get <- function(probes_prefix, grp)
{
  if(probes_prefix=="PROBES_CHR_")
  {
    probes <- semseeker::PROBES_CHR_CHR
    # message("DEBUG: loaded probes: PROBES_CHR_CHR")
  }
  else if(probes_prefix=="PROBES")
  {
    probes <- semseeker::PROBES
    # message("DEBUG: loaded probes: PROBES_CHR_CHR")
  }
  else
  {
    probes_name <- paste0(probes_prefix, grp,sep="")
    probes <- get(probes_name, envir = asNamespace("semseeker"))
    # message("DEBUG: loaded probes:", probes_name)
  }

  return(probes)
}
