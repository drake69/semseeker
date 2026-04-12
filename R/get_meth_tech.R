#' Detect the Illumina methylation array technology from a signal matrix
#'
#' Identifies whether a methylation signal matrix originates from a 27k, 450k,
#' or EPIC 850k array (or WGBS data).  Detection is performed in two steps:
#' \enumerate{
#'   \item \strong{Probe-ID matching}: probe identifiers are looked up against
#'     the bundled \code{PROBES} annotation table (columns \code{k27},
#'     \code{k450}, \code{k850}).  The technology with the most matching probes
#'     wins.  This approach is robust to subsetting (e.g. only 20 k probes out
#'     of 866 k).
#'   \item \strong{Row-count fallback}: used when \code{PROBES} is not
#'     available or yields no matches (e.g. WGBS data with non-cg identifiers).
#' }
#' The detected technology and beta/M-value flag are stored in the session
#' environment.
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

  ssEnv    <- get_session_info()
  n_probes <- nrow(signal_data)

  # Log informational row-count hints
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

  tech <- ""
  msg  <- ""

  # ---- Step 1: probe-ID matching against PROBES lookup table ----
  # PROBES is still bundled (transition period); columns are lowercase k27/k450/k850.
  # Once the large .rda files are removed this block will be skipped gracefully.
  if (tech == "" &&
      exists("PROBES", envir = asNamespace("SEMseeker"), inherits = FALSE)) {

    probe_map <- get("PROBES", envir = asNamespace("SEMseeker"))

    # Normalise column names: accept both k450/k850/k27 and K450/K850/K27
    col_map <- c(K27 = "k27", K450 = "k450", K850 = "k850")
    avail   <- intersect(c("k27", "k450", "k850", "K27", "K450", "K850"),
                         colnames(probe_map))
    if (length(avail) >= 3) {
      matched  <- probe_map[probe_map$PROBE %in% probe_ids, , drop = FALSE]
      if (nrow(matched) > 0) {
        # Count matches per technology (handle both naming conventions)
        counts <- setNames(integer(3), c("K27", "K450", "K850"))
        for (tname in c("K27", "K450", "K850")) {
          col_lower <- tolower(tname)
          col_upper <- tname
          col <- if (col_lower %in% colnames(matched)) col_lower else col_upper
          if (col %in% colnames(matched))
            counts[[tname]] <- sum(matched[[col]], na.rm = TRUE)
        }
        if (max(counts) > 0) {
          tech <- names(which.max(counts))
          msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                        "technology identified as", tech,
                        "from probe-ID matching (", max(counts), "probes matched).")
        }
      }
    }
  }

  # ---- Step 2: row-count fallback ----
  if (tech == "") {
    if (n_probes > 866562) {
      tech <- "WGBS"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as WGBS (row count > 866k).")
    } else if (all(grepl("_", probe_ids[seq_len(min(1000, length(probe_ids)))]))) {
      tech <- "WGBS"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "probe IDs contain underscores — treating as WGBS dataset.")
    } else if (n_probes <= 30000) {
      tech <- "K27"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 27k array (row count heuristic).")
    } else if (n_probes <= 510000) {
      tech <- "K450"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as 450k array (row count heuristic).")
    } else {
      tech <- "K850"
      msg  <- paste("INFO:", format(Sys.time(), "%a %b %d %X %Y"),
                    "dataset identified as EPIC 850k array (row count heuristic).")
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

  # ---- Detect beta vs M-values ----
  signal_cols <- signal_data[
    seq_len(min(10000, nrow(signal_data))),
    !colnames(signal_data) %in% c("PROBE", "CHR", "K27", "K450", "K850",
                                  "k27", "k450", "k850"),
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
