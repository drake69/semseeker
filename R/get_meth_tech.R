#' Detect the Illumina methylation array technology from a signal matrix
#'
#' Identifies whether a methylation signal matrix originates from a 27k, 450k,
#' or EPIC 850k array (or WGBS data) by matching probe row-counts and probe ID
#' patterns against known technology thresholds.  The detected technology and
#' beta/M-value flag are stored in the session environment.
#'
#' @param signal_data A numeric matrix or \code{data.frame} with CpG probes as
#'   rows and samples as columns.  Row names must be probe identifiers (e.g.
#'   \code{cg00000029}) unless a \code{PROBE} column is present.
#'
#' @return The updated session environment (\code{ssEnv}) invisibly.  The
#'   caller can retrieve the detected technology via
#'   \code{get_session_info()$tech}.
#'
get_meth_tech <- function(signal_data) {

  ssEnv <- get_session_info()

  n_probes <- nrow(signal_data)

  # Informational row-count heuristics
  if (n_probes == 485512)
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches 450k dataset.")
  if (n_probes == 27578)
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches 27k dataset.")
  if (n_probes == 866562)
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count matches EPIC 850k dataset.")
  if (n_probes > 866562)
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "probe count exceeds EPIC 850k — treating as WGBS dataset.")

  # Resolve probe identifiers
  probe_ids <- if ("PROBE" %in% colnames(signal_data))
    signal_data$PROBE
  else
    rownames(signal_data)

  # Detect technology from probe count thresholds
  tech <- ""
  msg  <- ""

  if (n_probes > 866562) {
    tech <- "WGBS"
    msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                  "dataset identified as WGBS.")
  } else {
    # Match probe IDs against the bundled DMR annotation to infer technology;
    # fall back to row-count ranges when dmr_annotation is not informative.
    n_cg <- sum(grepl("^cg", probe_ids[seq_len(min(1000, length(probe_ids)))]))

    if (n_probes <= 27578 || (n_cg > 0 && n_probes <= 30000)) {
      tech <- "K27"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 27k array.")
    } else if (n_probes <= 500000 || (n_cg > 0 && n_probes <= 510000)) {
      tech <- "K450"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 450k array.")
    } else if (n_probes <= 870000) {
      tech <- "K850"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as EPIC 850k array.")
    }

    # Fallback: if probe IDs contain underscores it is likely WGBS
    if (tech == "" &&
        all(grepl("_", probe_ids[seq_len(min(1000, length(probe_ids)))]))) {
      tech <- "WGBS"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "probe IDs contain underscores — treating as WGBS dataset.")
    }
  }

  if (tech == "") {
    msg <- paste("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),
                 "could not determine array technology.")
    log_event(msg)
    stop(msg)
  }

  log_event(msg)
  ssEnv$tech <- tech

  # Determine whether values are beta (0-1) or M-values (unbounded)
  signal_cols <- signal_data[
    seq_len(min(10000, nrow(signal_data))),
    !colnames(signal_data) %in% c("PROBE", "CHR", "K27", "K450", "K850"),
    drop = FALSE
  ]
  max_data   <- max(abs(c(max(signal_cols, na.rm = TRUE),
                          min(signal_cols, na.rm = TRUE))))
  ssEnv$beta <- max_data <= 1

  if (ssEnv$beta)
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "values are beta (0-1).")
  else
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
              "values appear to be M-values.")

  ssEnv$probes_count <- n_probes
  update_session_info(ssEnv)

  log_event("JOURNAL: array technology set to:", ssEnv$tech)
  return(ssEnv)
}
