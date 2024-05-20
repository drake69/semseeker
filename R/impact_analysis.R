#' @export
gene_impact_analysis <- function(inference_details, adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas,
  study, significance,statistic_parameter, sql_condition="",path_dbs, phenolyzer_folder_bin,disease,
  phenolyser=F, webgestalt=F, pathfindr=F,
  result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{
  start_fresh <- FALSE
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)

  arguments <- list(...)

  # check adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas have the same length
  if (length(adjust_per_area_s) != length(adjust_globally_s) | length(adjust_per_area_s) != length(pvalue_columns) | length(adjust_per_area_s) != length(adjustment_methods) | length(adjust_per_area_s) != length(alphas))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y")  ," adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas must have the same length !")
    stop()
  }

  for (id in 1:nrow(inference_details))
  {
    inference_detail <- inference_details[id,]
    for (i in 1:length(pvalue_columns)){

      pvalue_column <- pvalue_columns[i]
      alpha <- alphas[i]
      ssEnv$alpha <- alpha
      update_session_info(ssEnv)

      adjustment_method <- adjustment_methods[i]
      adjust_per_area <- adjust_per_area_s[i]
      adjust_globally <- adjust_globally_s[i]

      #
      if (webgestalt)
        pathway_WebGestalt(
          study = study,
          types=c("BP","MF"),
          enrich_methods = c("ORA"),
          adjust_per_area = adjust_per_area,
          adjust_globally = adjust_globally,
          adjustment_method = adjustment_method,
          pvalue_column = pvalue_column,
          inference_details = inference_detail,
          sql_condition = sql_condition,
          significance = significance
        )

      if(pathfindr)
        pathway_pathfindR(
          study =  study,
          path_db = path_dbs,
          iterations = 10,
          statistic_parameter=statistic_parameter,
          adjust_per_area = adjust_per_area,
          adjust_globally = adjust_globally,
          adjustment_method = adjustment_method,
          pvalue_column = pvalue_column,
          inference_details = inference_detail,
          sql_condition = sql_condition,
          significance = significance
        )

      if(phenolyser)
        phenotype_phenolyzer(
          study =  study,
          disease = disease,
          phenolyzer_folder_bin = phenolyzer_folder_bin,
          minimum_score = 0.5,
          statistic_parameter = statistic_parameter,
          adjust_per_area = adjust_per_area,
          adjust_globally = adjust_globally,
          adjustment_method = adjustment_method,
          pvalue_column = pvalue_column,
          inference_details = inference_detail,
          sql_condition = sql_condition,
          significance = significance
        )

    }
  }
}
