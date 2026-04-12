#' Build the probe annotation table from Bioconductor packages
#'
#' Constructs a \code{data.frame} equivalent to the legacy \code{pp_tot} object
#' by reading coordinates and feature annotations from the appropriate Illumina
#' array annotation package (installed via \code{minfi}).  The result is cached
#' inside the session environment (\code{ssEnv$probe_annotation}) so the
#' Bioconductor package is queried only once per R session.
#'
#' Column mapping from Bioconductor to SEMseeker:
#' \tabular{ll}{
#'   \code{PROBE}            \tab rownames of the annotation object         \cr
#'   \code{CHR}              \tab \code{chr} column (strip "chr" prefix)    \cr
#'   \code{START / END}      \tab \code{pos} column (single-base probes)    \cr
#'   \code{GENE_*}           \tab parsed from \code{UCSC_RefGene_Group}     \cr
#'   \code{ISLAND_*}         \tab recoded from \code{Relation_to_Island}    \cr
#'   \code{CHR_CYTOBAND}     \tab \code{Methyl450_Loci} or similar          \cr
#'   \code{DMR_WHOLE/DMR_DMR}\tab from bundled \code{dmr_annotation} table  \cr
#' }
#'
#' @param tech Character scalar: one of \code{"K27"}, \code{"K450"}, or
#'   \code{"K850"}.
#' @param force Logical: if \code{TRUE}, rebuild the cache even if it already
#'   exists.
#'
#' @return A \code{data.frame} with one row per probe and columns
#'   \code{PROBE}, \code{CHR}, \code{START}, \code{END}, plus feature columns
#'   for gene body regions, CpG islands, cytobands, and DMR annotations.
#'
#' @importFrom utils data
probe_annotation_build <- function(tech, force = FALSE) {

  ssEnv <- get_session_info()

  # Return cached version unless forced rebuild
  if (!force && !is.null(ssEnv$probe_annotation) &&
      !is.null(ssEnv$probe_annotation_tech) &&
      ssEnv$probe_annotation_tech == tech) {
    return(ssEnv$probe_annotation)
  }

  pkg <- switch(tech,
    K850 = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    K450 = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    K27  = "IlluminaHumanMethylation27kanno.ilmn12.hg19",
    stop("Unknown array technology: ", tech,
         ". Must be one of K27, K450, K850.")
  )

  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Annotation package '", pkg, "' is required but not installed.\n",
         "Install it with: BiocManager::install('", pkg, "')")

  if (!requireNamespace("minfi", quietly = TRUE))
    stop("Package 'minfi' is required. Install with: BiocManager::install('minfi')")

  anno_obj <- get(pkg, envir = asNamespace(pkg))
  anno     <- minfi::getAnnotation(anno_obj)
  anno_df  <- as.data.frame(anno, stringsAsFactors = FALSE)

  # ---- PROBE, CHR, START, END ----
  anno_df$PROBE <- rownames(anno_df)
  anno_df$CHR   <- sub("^chr", "", anno_df$chr)   # "chr1" -> "1"
  anno_df$START <- anno_df$pos
  anno_df$END   <- anno_df$pos

  # ---- Technology flag column ----
  anno_df[[tech]] <- TRUE

  # ---- GENE columns ----
  # UCSC_RefGene_Group: semicolon-separated list of region types per transcript
  # UCSC_RefGene_Name:  semicolon-separated list of gene symbols
  gene_groups <- strsplit(as.character(anno_df$UCSC_RefGene_Group), ";")
  gene_names  <- strsplit(as.character(anno_df$UCSC_RefGene_Name),  ";")

  gene_region_cols <- c(
    GENE_BODY     = "Body",
    GENE_TSS200   = "TSS200",
    GENE_TSS1500  = "TSS1500",
    GENE_1STEXON  = "1stExon",
    GENE_5UTR     = "5'UTR",
    GENE_3UTR     = "3'UTR",
    GENE_EXONBND  = "ExonBnd"
  )

  for (col in names(gene_region_cols)) {
    region <- gene_region_cols[[col]]
    anno_df[[col]] <- mapply(function(groups, genes) {
      hits <- genes[groups == region]
      if (length(hits) == 0) return(NA_character_)
      paste(unique(hits[hits != ""]), collapse = ";")
    }, gene_groups, gene_names, SIMPLIFY = TRUE)
  }

  # GENE_WHOLE: any gene associated with the probe
  anno_df$GENE_WHOLE <- mapply(function(genes) {
    hits <- unique(genes[genes != "" & !is.na(genes)])
    if (length(hits) == 0) return(NA_character_)
    paste(hits, collapse = ";")
  }, gene_names, SIMPLIFY = TRUE)

  # ---- ISLAND columns ----
  # Relation_to_Island values: Island, N_Shore, S_Shore, N_Shelf, S_Shelf, OpenSea
  island_rel <- anno_df$Relation_to_Island

  anno_df$ISLAND_WHOLE    <- ifelse(island_rel == "Island",
                                    anno_df$Islands_Name, NA_character_)
  anno_df$ISLAND_N_SHORE  <- ifelse(island_rel == "N_Shore",
                                    anno_df$Islands_Name, NA_character_)
  anno_df$ISLAND_S_SHORE  <- ifelse(island_rel == "S_Shore",
                                    anno_df$Islands_Name, NA_character_)
  anno_df$ISLAND_N_SHELF  <- ifelse(island_rel == "N_Shelf",
                                    anno_df$Islands_Name, NA_character_)
  anno_df$ISLAND_S_SHELF  <- ifelse(island_rel == "S_Shelf",
                                    anno_df$Islands_Name, NA_character_)

  # ---- CHR_CYTOBAND ----
  if ("UCSC_CpG_Islands_Name" %in% colnames(anno_df))
    anno_df$CHR_CYTOBAND <- anno_df$UCSC_CpG_Islands_Name
  else if ("Methyl450_Loci" %in% colnames(anno_df))
    anno_df$CHR_CYTOBAND <- anno_df$Methyl450_Loci
  else
    anno_df$CHR_CYTOBAND <- NA_character_

  # ---- DMR columns — from bundled dmr_annotation ----
  dmr_annotation <- SEMseeker::dmr_annotation
  anno_df <- merge(anno_df, dmr_annotation, by = "PROBE", all.x = TRUE)

  # ---- Select and return only SEMseeker columns ----
  keep_cols <- c(
    "PROBE", "CHR", "START", "END", tech,
    "GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
    "GENE_5UTR", "GENE_3UTR", "GENE_EXONBND", "GENE_WHOLE",
    "ISLAND_WHOLE", "ISLAND_N_SHORE", "ISLAND_S_SHORE",
    "ISLAND_N_SHELF", "ISLAND_S_SHELF",
    "CHR_CYTOBAND", "DMR_WHOLE", "DMR_DMR"
  )
  keep_cols <- intersect(keep_cols, colnames(anno_df))
  anno_df   <- anno_df[, keep_cols, drop = FALSE]

  # Cache in session environment
  ssEnv$probe_annotation      <- anno_df
  ssEnv$probe_annotation_tech <- tech
  update_session_info(ssEnv)

  log_event("INFO: probe annotation built from Bioconductor package '", pkg,
            "' (", nrow(anno_df), " probes)")

  return(anno_df)
}
