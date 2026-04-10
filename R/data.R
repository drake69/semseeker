#' PROBES_CHR_CHR
#'
#' Probe-to-chromosome mapping table for Illumina methylation arrays.
#' Maps each CpG probe to its chromosome-level annotation.
#'
#' @format A data frame with columns for probe ID, chromosome, start/end
#'   position, and chromosome-level annotation features.
"PROBES_CHR_CHR"


#' PROBES
#'
#' Subset of probe annotations for Illumina methylation arrays (450k/EPIC).
#' Contains genomic coordinates and feature annotations for each CpG probe.
#'
#' @format A data frame with columns including \code{PROBE} (probe identifier),
#'   \code{CHR} (chromosome), \code{START} (genomic position), \code{END},
#'   \code{K27}, \code{K450}, \code{K850} (array technology flags), and
#'   additional gene/island annotation columns.
"PROBES"


#' pp_tot
#'
#' Full probe annotation table for Illumina methylation arrays (27k, 450k, EPIC 850k).
#' Used internally by \code{probe_features_get()} to look up genomic coordinates
#' and feature annotations for each CpG probe.
#'
#' @format A data frame with approximately 800,000 rows (one per CpG probe) and
#'   columns including \code{PROBE}, \code{CHR}, \code{START}, \code{END},
#'   \code{K27}, \code{K450}, \code{K850} (array technology flags), gene body
#'   feature columns (\code{GENE_BODY}, \code{GENE_TSS200}, etc.), and CpG
#'   island context columns (\code{ISLAND_N_SHORE}, \code{ISLAND_S_SHORE}, etc.).
"pp_tot"


#' metrics_properties
#'
#' Metadata table describing the statistical properties of each SEM metric.
#' Used internally to determine ranking direction and scaling behaviour during
#' association analysis and pathway enrichment.
#'
#' @format A data frame with one row per metric and columns including
#'   \code{Metric} (metric name), \code{Higher_the_Better} (logical: whether
#'   higher values indicate stronger signal), and \code{Affected_by_Scaling}
#'   (logical: whether the metric is affected by data transformation).
"metrics_properties"


#' ssEnv
#'
#' Internal session environment object persisted between SEMseeker analysis
#' steps. Stores runtime parameters (result folder paths, technology flag,
#' alpha threshold, etc.) set by \code{\link{init_env}} and retrieved by
#' \code{get_session_info()}.
#'
#' @format An \code{environment} containing named slots for session-level
#'   analysis parameters. Users should not modify this object directly;
#'   use \code{\link{init_env}} and \code{set_env_variable()} instead.
"ssEnv"
