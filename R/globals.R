# Suppress R CMD check NOTEs for variables used in data-masking contexts
# (dplyr, ggplot2, data.table) and for package-level globals.

utils::globalVariables(c(
  # Column names used in dplyr/ggplot2/data.table masking
  "AREA", "DEPTH", "ID", "MARKER", "P_Value", "Q", "SAMPLE_GROUP",
  "SOURCE", "X1", "X2", "covariates", "data", "db", "k", "log_fdr",
  "pathway_report", "results", "study", "test.data",
  "value", "variable", "x",
  # Package session environment used across functions
  "ssEnv"
))
