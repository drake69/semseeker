#' Build the probe annotation table from Bioconductor packages
#'
#' Constructs a \code{data.frame} equivalent to the legacy \code{pp_tot} object
#' by reading coordinates and feature annotations directly from the S4 data
#' slots of the installed Illumina annotation package.  \pkg{minfi} is \emph{not}
#' required at runtime; the annotation packages are accessed without it.
#' The result is cached inside the session environment
#' (\code{ssEnv$probe_annotation}) so the package is parsed only once per
#' R session.
#'
#' Column mapping from Bioconductor to SEMseeker:
#' \tabular{ll}{
#'   \code{PROBE}             \tab rownames of the annotation Locations table  \cr
#'   \code{CHR}               \tab \code{chr} column (strip "chr" prefix)      \cr
#'   \code{START / END}       \tab \code{pos} column (single-base CpG probes)  \cr
#'   \code{GENE_*}            \tab parsed from \code{UCSC_RefGene_Group}       \cr
#'   \code{ISLAND_*}          \tab recoded from \code{Relation_to_Island}      \cr
#'   \code{CHR_CYTOBAND}      \tab \code{UCSC_CpG_Islands_Name} or similar     \cr
#'   \code{DMR_WHOLE/DMR_DMR} \tab from bundled \code{dmr_annotation} table    \cr
#' }
#'
#' @param tech Character scalar: one of \code{"K27"}, \code{"K450"}, or
#'   \code{"K850"}.
#' @param force Logical: if \code{TRUE}, rebuild the cache even when a cached
#'   version already exists.
#'
#' @return A \code{data.frame} with one row per probe and columns
#'   \code{PROBE}, \code{CHR}, \code{START}, \code{END}, the technology flag,
#'   gene body region columns, CpG island context columns, cytoband, and DMR
#'   annotations.
#'
probe_annotation_build <- function(tech, force = FALSE) {

  ssEnv <- get_session_info()

  # Return cached version unless forced
  if (!force &&
      !is.null(ssEnv$probe_annotation) &&
      identical(ssEnv$probe_annotation_tech, tech)) {
    return(ssEnv$probe_annotation)
  }

  pkg <- .ANNO_PKGS[[tech]]
  if (is.null(pkg))
    stop("Unknown technology: '", tech, "'. Must be one of: ",
         paste(names(.ANNO_PKGS), collapse = ", "))

  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Annotation package '", pkg, "' is required but not installed.\n",
         "Install with: BiocManager::install('", pkg, "')")

  log_event("INFO: building probe annotation from '", pkg, "'...")

  # ---- Read raw data from S4 slots (no minfi required) ----
  anno_df <- .anno_pkg_to_df(pkg)

  # ---- PROBE, CHR, START, END ----
  anno_df$PROBE <- rownames(anno_df)
  anno_df$CHR   <- sub("^chr", "", as.character(anno_df$chr))
  anno_df$START <- as.integer(anno_df$pos)
  anno_df$END   <- as.integer(anno_df$pos)

  # ---- Technology flag ----
  anno_df[[tech]] <- TRUE

  # ---- GENE columns ----
  gene_groups <- strsplit(as.character(anno_df$UCSC_RefGene_Group), ";")
  gene_names  <- strsplit(as.character(anno_df$UCSC_RefGene_Name),  ";")

  gene_region_map <- c(
    GENE_BODY    = "Body",
    GENE_TSS200  = "TSS200",
    GENE_TSS1500 = "TSS1500",
    GENE_1STEXON = "1stExon",
    GENE_5UTR    = "5'UTR",
    GENE_3UTR    = "3'UTR",
    GENE_EXONBND = "ExonBnd"
  )
  for (col in names(gene_region_map)) {
    region <- gene_region_map[[col]]
    anno_df[[col]] <- mapply(function(groups, genes) {
      hits <- unique(genes[groups == region & genes != ""])
      if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
    }, gene_groups, gene_names, SIMPLIFY = TRUE)
  }

  anno_df$GENE_WHOLE <- mapply(function(genes) {
    hits <- unique(genes[genes != "" & !is.na(genes)])
    if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
  }, gene_names, SIMPLIFY = TRUE)

  # ---- ISLAND columns ----
  island_rel <- as.character(anno_df$Relation_to_Island)
  island_name <- as.character(anno_df$Islands_Name)

  anno_df$ISLAND_WHOLE   <- ifelse(island_rel == "Island",  island_name, NA_character_)
  anno_df$ISLAND_N_SHORE <- ifelse(island_rel == "N_Shore", island_name, NA_character_)
  anno_df$ISLAND_S_SHORE <- ifelse(island_rel == "S_Shore", island_name, NA_character_)
  anno_df$ISLAND_N_SHELF <- ifelse(island_rel == "N_Shelf", island_name, NA_character_)
  anno_df$ISLAND_S_SHELF <- ifelse(island_rel == "S_Shelf", island_name, NA_character_)

  # ---- CHR_CYTOBAND ----
  cytoband_col <- intersect(
    c("UCSC_CpG_Islands_Name", "Methyl450_Loci", "Methyl27_Loci"),
    colnames(anno_df)
  )
  anno_df$CHR_CYTOBAND <- if (length(cytoband_col) > 0L)
    as.character(anno_df[[cytoband_col[[1L]]]])
  else
    NA_character_

  # ---- DMR columns from bundled dmr_annotation ----
  dmr <- SEMseeker::dmr_annotation
  anno_df <- merge(anno_df, dmr, by = "PROBE", all.x = TRUE)

  # ---- Select final columns ----
  keep <- c(
    "PROBE", "CHR", "START", "END", tech,
    "GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
    "GENE_5UTR", "GENE_3UTR", "GENE_EXONBND", "GENE_WHOLE",
    "ISLAND_WHOLE", "ISLAND_N_SHORE", "ISLAND_S_SHORE",
    "ISLAND_N_SHELF", "ISLAND_S_SHELF",
    "CHR_CYTOBAND", "DMR_WHOLE", "DMR_DMR"
  )
  anno_df <- anno_df[, intersect(keep, colnames(anno_df)), drop = FALSE]

  # ---- Cache and return ----
  ssEnv$probe_annotation      <- anno_df
  ssEnv$probe_annotation_tech <- tech
  update_session_info(ssEnv)

  log_event("INFO: probe annotation built — ", nrow(anno_df), " probes, tech = ", tech)
  return(anno_df)
}
