coverage_analysis <- function(methylation_data)
{
  # probes <- PROBES_CHR_CHR
  # grp <- c("CHR")
  # probes_prefix <- "PROBES_CHR_"

  keys <- rbind(
    envir$keys_gene_probes,
    envir$keys_dmr_probes,
    envir$keys_island_probes,
    # envir$keys_probe_probes,
    envir$keys_chr_probes
  )

  for ( k in 1:nrow(keys))
  {
    # k <- 16
    SUBGROUP <- as.character(keys[k,"subgroups"])
    probes_prefix <- as.character(keys[k,"prefix"])
    grp <-  as.character(keys[k,"maingrouplable"])
    # message(grp,"\n")

    probes <- probes_get(probes_prefix, SUBGROUP)
    probe_filtered <- probes[ probes$PROBE %in% rownames(methylation_data),]
    covered_count <- aggregate(probe_filtered$PROBE,list(probe_filtered[,grp]), FUN=length)

    colnames(covered_count) <- c(as.character(grp),"COUNT_COVERED")
    total_count <- aggregate(probes$PROBE, list(probes[, grp]), FUN=length)
    colnames(total_count) <- c(as.character(grp),"COUNT_TOTAL")

    coverage <- merge(total_count, covered_count, by=grp, all.x = TRUE)
    coverage <- coverage[ !is.na(coverage[,grp]),]
    coverage <- coverage[ coverage[,grp]!="",]
    coverage$PERC <- plyr::round_any(100 * coverage$COUNT_COVERED / coverage$COUNT_TOTAL, 10)

    cov_stat <- aggregate(coverage$PERC, list(coverage$PERC), FUN=length)
    colnames(cov_stat) <- c("COV_PERC","COUNT")
    cov_stat$GROUP <- grp
    cov_stat$SUBGROUP <- SUBGROUP
    if(exists("cov_result"))
      cov_result <- rbind(cov_result, cov_stat)
    else
      cov_result <- cov_stat
  }

}
