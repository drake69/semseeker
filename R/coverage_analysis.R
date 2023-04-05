coverage_analysis <- function(methylation_data, envir)
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
    coverage$PERC <- plyr::round_any(100 * coverage$COUNT_COVERED / coverage$COUNT_TOTAL, 5)

    cov_stat <- aggregate(coverage$PERC, list(coverage$PERC), FUN=length)
    colnames(cov_stat) <- c("COV_PERC","COUNT")
    cov_stat$GROUP <- grp
    cov_stat$SUBGROUP <- SUBGROUP
    if(exists("cov_result"))
      cov_result <- rbind(cov_result, cov_stat)
    else
      cov_result <- cov_stat

    colnames(total_count) <- c("GROUP","COUNT_TOTAL")
    total_count$AREA <- grp
    total_count$SUBGROUP <- SUBGROUP

    if(exists("tot_result"))
      tot_result <- rbind(total_count, tot_result)
    else
      tot_result <- total_count
  }

  cov_result <- reshape2::dcast(data = cov_result, GROUP + SUBGROUP ~ COV_PERC, value.var = "COUNT", sum)
  tot_result <- subset(tot_result, SUBGROUP!="CHR" & SUBGROUP != "Whole")
  # tot_result <- aggregate(tot_result$SUBGROUP, list(tot_result$COUNT_TOTAL), FUN=length)
  tot_result <- reshape2::dcast(data = tot_result, AREA + SUBGROUP ~ COUNT_TOTAL, value.var = "COUNT_TOTAL", length)

  chartFolder <- dir_check_and_create(envir$result_folderChart,"COVERAGE")
  tt <- as.data.frame(cov_result[,3:ncol(cov_result)])
  rownames(tt) <- paste(cov_result$GROUP,cov_result$SUBGROUP, sep=" ")
  colors <- grDevices::hsv(0.560, seq(0,1,length.out = 20) , 1)
  if (!plyr::empty(tt))
    if(nrow(tt) > 2 & ncol(tt) > 2)
    {
      filename = paste0( chartFolder,"/","COVERAGE_ANALYSIS.png",sep="")
      grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = 144)
      stats::heatmap(
        x =  as.matrix(tt),
        col = colors,
        scale = "column",
        Colv = NA,
        Rowv = NA,
        # RowSideColors = as.vector( as.character(cov_result$COV_PERC)),
        margins = c(15, 15),
        main = "Coverage analysis of probes for each genomic area"
      )
      graphics::legend(x="right", legend=c("min", "med", "max"), fill=colors[c(1,10,20)])
      grDevices::dev.off()
    }

  return(cov_result)
}
