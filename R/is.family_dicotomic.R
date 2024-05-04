is.family_dicotomic <- function(family_test) {

  if(family_test=="gaussian" | family_test=="spearman" | family_test=="pearson" |
      family_test=="kendall" | grepl("quantreg-permutation", family_test) | grepl("quantreg", family_test)
    | family_test=="poisson"  | grepl("mean-permutation", family_test) |
      grepl("polynomial", family_test) | grepl("exp", family_test)  | grepl("log", family_test))
    return(FALSE)

  if(family_test=="wilcoxon" | family_test=="t.test" | family_test =="jsd" | family_test=="chisq.test"
    | family_test=="fisher.test" | family_test=="kruskal.test" | family_test=="binomial")
    return(TRUE)

  log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " family_test not recognized ")
  stop()
}
