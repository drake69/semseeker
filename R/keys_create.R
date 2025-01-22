keys_create <- function(ssEnv, arguments)
{

  # remove BOTH and BOTHSUM
  #
  ssEnv$keys_sample_groups <-  data.frame("SAMPLE_GROUP"=c("Reference","Control","Case"))
  ##### MANAGE MARKERS AND FIGURES #####
  # set default values
  keys_figures_default_discrete <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH","BOTHSUM")) #
  keys_markers_default_discrete <-  data.frame("MARKER"=c("MUTATIONS","LESIONS","DELTAQ","DELTARQ","DELTAP","DELTARP"))


  keys_markers_default_discrete$Q <- c(1,NA,ssEnv$DELTAQ_Q,ssEnv$DELTARQ_Q,ssEnv$DELTAP_B,ssEnv$DELTARP_B)

  keys_markers_default_discrete$SUFFIX <-  rep("",6)
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_discrete$EXT <-  c("bed","bed","bed","bed","bed","bed")
  keys_markers_figures_default_discrete <- merge(keys_markers_default_discrete, keys_figures_default_discrete, by = NULL)
  rm(keys_figures_default_discrete, keys_markers_default_discrete)


  keys_figures_default_continuos <-  data.frame("FIGURE"=c("HYPO", "HYPER", "BOTH","BOTHSUM")) # , "BOTH","BOTHSUM"
  keys_markers_default_continuos <-  data.frame("MARKER"=c("DELTAS","DELTAR"))
  keys_markers_default_continuos$SUFFIX <-  c("","")
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_continuos$EXT <-  c("bedgraph","bedgraph")
  keys_markers_figures_default_continuos <- merge(keys_markers_default_continuos, keys_figures_default_continuos, by = NULL)
  rm(keys_figures_default_continuos, keys_markers_default_continuos)

  keys_figures_default_continuos_mono_figure <-  data.frame("FIGURE"=c("MEAN"))
  keys_markers_default_continuos_mono_figure <-  data.frame("MARKER"=c("SIGNAL"))
  keys_markers_default_continuos_mono_figure$SUFFIX <-  c("")
  # keys_markers_default$SUFFIX <-  c("","","","", ssEnv$epiquantile ,"",ssEnv$epiquantile)
  keys_markers_default_continuos_mono_figure$EXT <-  c("bedgraph")
  keys_markers_figures_default_continuos_mono_figure <- merge(keys_markers_default_continuos_mono_figure, keys_figures_default_continuos_mono_figure, by = NULL)
  rm(keys_figures_default_continuos_mono_figure, keys_markers_default_continuos_mono_figure)

  keys_markers_figures_default <- plyr::rbind.fill(keys_markers_figures_default_discrete,
    keys_markers_figures_default_continuos,
    keys_markers_figures_default_continuos_mono_figure)

  ssEnv$keys_markers_figures_default <- keys_markers_figures_default

  figures <- arguments[["figures"]]
  arguments[["figures"]] <- NULL
  if (!is.null(figures))
    keys_markers_figures_default <- keys_markers_figures_default[keys_markers_figures_default$FIGURE %in% figures,]

  markers <- arguments[["markers"]]
  arguments[["markers"]] <- NULL
  if (!is.null(markers))
    keys_markers_figures_default <- keys_markers_figures_default[keys_markers_figures_default$MARKER %in% markers,]

  ########################################################################################################################
  # manage AREAS and SUBAREAS #####

  keys_areas_default <- data.frame("AREA"=c("GENE","ISLAND","DMR","CHR","PROBE"))

  keys_gene_subareas_default <- data.frame("AREA"="GENE", "SUBAREA"=c("BODY","TSS1500","TSS200","1STEXON","3UTR","5UTR","EXONBND","WHOLE"))
  keys_island_subareas_default <- data.frame("AREA"="ISLAND","SUBAREA"=c("N_SHORE","S_SHORE","N_SHELF","S_SHELF", "WHOLE"))
  keys_dmr_subareas_default <- data.frame("AREA"="DMR","SUBAREA"=c("WHOLE","DMR"))
  keys_chr_subareas_default <- data.frame("AREA"="CHR","SUBAREA"=c("WHOLE","CYTOBAND"))
  keys_probe_subareas_default <- data.frame("AREA"="PROBE","SUBAREA"=c("WHOLE"))

  keys_areas_subareas_default <- rbind(keys_gene_subareas_default, keys_island_subareas_default, keys_dmr_subareas_default, keys_chr_subareas_default, keys_probe_subareas_default)

  areas <- arguments[["areas"]]
  arguments[["areas"]] <- NULL
  if (!is.null(areas))
    keys_areas_subareas_default <- keys_areas_subareas_default[keys_areas_subareas_default$AREA %in% areas,]

  subareas <- arguments[["subareas"]]
  arguments[["subareas"]] <- NULL
  if (!is.null(subareas))
    keys_areas_subareas_default <- keys_areas_subareas_default[keys_areas_subareas_default$SUBAREA %in% subareas,]

  ########################################################################################################################
  # combine AREAS, SUBAREAS, MARKERS and FIGURES #####

  keys_areas_subareas_markers_figures <- unique(merge(keys_areas_subareas_default, keys_markers_figures_default, by = NULL))

  ########################################################################################################################

  # assign to ssEnv
  ssEnv$keys_areas_subareas <- keys_areas_subareas_default
  ssEnv$keys_markers_figures <- keys_markers_figures_default
  ssEnv$keys_areas_subareas_markers_figures <- keys_areas_subareas_markers_figures

  combine_not_empty <- function(x)
  {
    # paste0(x[x!=""], collapse = "_")
    # remove spaces
    gsub(" ","",paste0(x[x!=""], collapse = "_"))
  }

  ssEnv$keys_areas_subareas_markers_figures$COMBINED <- apply(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","FIGURE","AREA","SUBAREA")], 1, combine_not_empty )
  ssEnv$keys_areas_subareas$COMBINED <- apply(ssEnv$keys_areas_subareas[,c("AREA","SUBAREA")], 1, combine_not_empty )
  ssEnv$keys_markers_figures$COMBINED <- apply(ssEnv$keys_markers_figures[,c("MARKER","FIGURE")], 1, combine_not_empty )

  ########################################################################################################################
  # prepare keys for pathway analysis
  keys <- ssEnv$keys_areas_subareas_markers_figures
  # select only gene
  keys <- keys[keys$AREA=="GENE",]
  keys$FIGURE <- as.character(keys$FIGURE)
  keys$SUBAREA <- as.character(keys$SUBAREA)
  # remove "COMBINED"
  keys <- keys[,c("AREA","SUBAREA","MARKER","FIGURE")]
  # where FIGURE is HYPO or HYPER for FIGURE set to HYPER_HYPO
  selector <- (keys$FIGURE=="HYPER" | keys$FIGURE=="HYPO")
  if (any(selector))
    keys[selector,"FIGURE"] <- "HYPER_HYPO"
  keys <- unique(keys)
  parts <- ssEnv$keys_areas_subareas[ssEnv$keys_areas_subareas[,"SUBAREA"]!="WHOLE","SUBAREA"]
  keys_parts <- keys[keys$SUBAREA %in% parts,]
  if (nrow(keys_parts)>0)
    keys_parts$SUBAREA <- "ALL_SUBAREAS"
  keys_whole <- keys[keys$SUBAREA=="WHOLE",]
  keys <- rbind(keys_parts,keys_whole)
  keys <- keys[!duplicated(keys),]
  # remove only HYPER or only HYPO
  keys <- keys[!(keys$FIGURE=="HYPER" | keys$FIGURE=="HYPO"),]
  keys$COMBINED <- paste(keys$AREA,keys$SUBAREA,keys$MARKER,keys$FIGURE, sep="_")
  # remove BOTH and BOTHSUM
  keys <- keys[!(keys$FIGURE=="BOTH" | keys$FIGURE=="BOTHSUM"),]
  ssEnv$keys_for_pathway <- keys

  ########################################################################################################################
  # prepare keys for enrichment analysis report format
  key_enrichment_format <- data.frame("type"="Pathway", "label"="STRINGdb","column_of_id"="term","column_of_description"="description", "column_of_pvalue"="fdr","column_of_enrichment"="fold_enrichment")
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="Phenolyzer_STRINGdb","column_of_id"="term","column_of_description"="description", "column_of_pvalue"="fdr","column_of_enrichment"="fold_enrichment"))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Pathway","label"="WebGestalt","column_of_id"="geneSet","column_of_description"="description","column_of_pvalue"="FDR","column_of_enrichment"="enrichmentRatio"))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Pathway","label"="Phenolyzer_WebGestalt","column_of_id"="geneSet","column_of_description"="description","column_of_pvalue"="FDR","column_of_enrichment"="enrichmentRatio"))
  key_enrichment_format <- rbind(key_enrichment_format, data.frame("type"="Phenotype","label"="phenolyzer","column_of_id"="Description","column_of_description"="Description","column_of_pvalue"="Score","column_of_enrichment"="Score"))
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="pathfindR","column_of_id"="ID","column_of_description"="Term_Description", "column_of_pvalue"="highest_p","column_of_enrichment"="Fold_Enrichment"))
  key_enrichment_format <- rbind(key_enrichment_format,data.frame("type"="Pathway", "label"="ctdR","column_of_id"="ChemicalID","column_of_description"="ChemicalName", "column_of_pvalue"="padj","column_of_enrichment"="foldEnrichment"))

  ssEnv$key_enrichment_format <- key_enrichment_format


  update_session_info(ssEnv)

  return(arguments)
}
