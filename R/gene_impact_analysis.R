#' @export
gene_impact_analysis <- function(inference_details, adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas,
  study, significance,statistic_parameter, path_dbs, phenolyzer_folder_bin,disease,
  phenolyser=FALSE, WebGestalt=FALSE, pathfindr=FALSE,STRINGdb=FALSE,Phenolyzer_STRINGdb=FALSE,Phenolyzer_WebGestalt=FALSE,ctdR=FALSE,
  result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{
  start_fresh <- FALSE
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = start_fresh, ...)

  arguments <- list(...)

  # check adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods,alphas have the same length
  if (length(adjust_per_area_s) != length(adjust_globally_s)
    | length(adjust_per_area_s) != length(pvalue_columns)
    | length(adjust_per_area_s) != length(adjustment_methods))
  {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y")  ," adjust_per_area_s, adjust_globally_s, pvalue_columns, adjustment_methods must have the same length !")
    stop()
  }


  inference_details <- subset(inference_details, depth_analysis ==3)
  for ( a in alphas)
  {
    for (id in 1:nrow(inference_details))
    {
      inference_detail <- inference_details[id,]
      inference_detail_prettified <- t(inference_detail)
      # Generate a plain text table using kable
      inference_detail_prettified <- knitr::kable(inference_detail_prettified, format = "simple",
        align = "l",    # Left align for all columns
        digits = 2,     # Number of digits for numeric columns
        row.names = TRUE) # Suppress row names

      inference_detail_prettified <- paste(inference_detail_prettified, collapse = "\n")
      log_event("JOURNAL: ##############################################################################################################")
      log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " \nStarting pathway for inference detail:\n", inference_detail_prettified)
      log_event("JOURNAL: Gene are selected by the sql_area_selection in the inference details above, also filtered using ", pvalue_columns, "columns. \nWith an alpha limit having value: ", a)

      for (i in seq_along(pvalue_columns)){

        pvalue_column <- pvalue_columns[i]
        ssEnv$alpha <- a
        update_session_info(ssEnv)

        adjustment_method <- adjustment_methods[i]
        adjust_per_area <- adjust_per_area_s[i]
        adjust_globally <- adjust_globally_s[i]

        if(ctdR)
          pathway_ctdR(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
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
            inference_detail = inference_detail,
            significance = significance
          )

        if(Phenolyzer_STRINGdb)
          pathway_Phenolyzer_STRINGdb(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance,
            disease = disease
          )


        if (Phenolyzer_WebGestalt)
          pathway_Phenolyzer_WebGestalt(
            study = study,
            types=c("BP","MF"),
            enrich_methods = c("ORA"),
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance,
            disease = disease
          )


        if(STRINGdb)
          pathway_STRINGdb(
            study =  study,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

        if(pathfindr)
          pathway_pathfindR(
            study =  study,
            path_db = path_dbs,
            iterations = 20,
            statistic_parameter=statistic_parameter,
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_details = inference_detail,
            significance = significance
          )


        if (WebGestalt)
          pathway_WebGestalt(
            study = study,
            types=c("BP","MF"),
            enrich_methods = c("ORA"),
            adjust_per_area = adjust_per_area,
            adjust_globally = adjust_globally,
            adjustment_method = adjustment_method,
            pvalue_column = pvalue_column,
            inference_detail = inference_detail,
            significance = significance
          )

      }
    }
  }

  close_env()
}
