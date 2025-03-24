execute_model <- function(family_test, tempDataFrame, sig.formula, burdenValue, independent_variable, transformation="", plot=FALSE,
  samples_sql_condition="", area="", subarea=""){

  samples_sql_condition = as.character(samples_sql_condition)
  transformation <- as.character(transformation)
  family_test <- as.character(family_test)
  # model_result <- data.frame()
  #
  # tryCatch({
  if( family_test=="binomial" | family_test=="poisson" | family_test=="gaussian")
    model_result <- glm_model(family_test, tempDataFrame, sig.formula, transformation, plot , samples_sql_condition=samples_sql_condition, area, subarea)

  if(family_test=="wilcoxon" | family_test=="t.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" | family_test=="jsd"
    | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test" | family_test=="bartlett.test")
    model_result <- test_model(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable, transformation, plot, samples_sql_condition, area, subarea)

  if (grepl("wilcoxon.paired",family_test) | grepl("t.test.paired",family_test ))
    model_result <- test_model_paired(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("mean-permutation",family_test))
    model_result <- mean_permutation(family_test, sig.formula, tempDataFrame, independent_variable,plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("quantile-permutation",family_test))
    model_result <- quantile_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("spearman-permutation",family_test))
    model_result <- spearman_permutation(family_test, sig.formula, tempDataFrame, independent_variable, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("polynomial",family_test))
    model_result <- polynomial_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("exp",family_test))
    model_result <- exp_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("log",family_test))
    model_result <- log_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("mediation-quantreg",family_test))
    model_result <- mediation_quantreg_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("mediation-ridge",family_test))
    model_result <- mediation_ridge_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if(grepl("mediation-linear",family_test))
    model_result <- mediation_linear_model(family_test, tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea)

  if (grepl("quantreg-permutation", family_test))
  {
    # Determine the null device for the current platform
    null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
    # Redirect output to the null device
    sink(null_device)
    model_result <- quantreg_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable, transformation, plot , samples_sql_condition=samples_sql_condition, area, subarea)
    sink()
  }
  if (grepl("quantreg", family_test) & !grepl("quantreg-permutation", family_test)  & !grepl("mediation-quantreg", family_test))
    model_result <- quantreg_model(family_test, sig.formula, tempDataFrame, independent_variable, transformation, plot,
      samples_sql_condition=samples_sql_condition, area, subarea)

  #   return (model_result)
  # }, error = function(e) {
  #   return (model_result)
  # })

  #
  model_result <- as.data.frame(model_result)
  return (model_result)

}
