execute_model <- function(family_test, tempDataFrame, sig.formula, burdenValue, independent_variable){

  if( family_test=="binomial" | family_test=="poisson" | family_test=="gaussian")
    model_result <- glm_model(family_test, tempDataFrame, sig.formula )

  if(family_test=="wilcoxon" | family_test=="t.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" | family_test=="jsd"
    | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test")
    model_result <- test_model(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable )

  if(grepl("mean-permutation",family_test))
    model_result <- mean_permutation(family_test, sig.formula, tempDataFrame, independent_variable)

  if(grepl("quantile-permutation",family_test))
    model_result <- quantile_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable)

  if(grepl("spearman-permutation",family_test))
    model_result <- spearman_permutation(family_test, sig.formula, tempDataFrame, independent_variable)

  if(grepl("polynomial",family_test))
    model_result <- polynomial_model(family_test, tempDataFrame, sig.formula)

  # Determine the null device for the current platform
  null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  # Redirect output to the null device
  sink(null_device)
  if (grepl("quantreg-permutation", family_test))
    model_result <- quantreg_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable)
  sink()
  # sink(file.path(ssEnv$session_folder,"session_output.log"), split = TRUE, append = TRUE)

  if (grepl("quantreg", family_test))
    model_result <- quantreg_model(family_test, sig.formula, tempDataFrame, independent_variable)

  return (model_result)
}
