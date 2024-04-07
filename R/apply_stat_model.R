#' Title
#'
#' @param tempDataFrame data frame to apply association
#' @param g_start index of starting data
#' @param family_test family of test to run
#' @param covariates vector of covariates
#' @param key key to identify file to elaborate
#' @param transformation transformation to apply to covariates, burden and independent variable
#' @param dototal do a total per area
#' @param session_folder where to save log file
#' @param independent_variable independent variable name
#' @param depth_analysis depth's analysis
#' @param ... extra parameters
#'
#' @importFrom doRNG %dorng%
#'
apply_stat_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation, dototal, session_folder,
                             independent_variable, depth_analysis=3, ...)
{
  ssEnv <- get_session_info()
  arguments <- list(...)

  browser()


  sig.formula <- as.formula("prova ~ Sample_Group")
  independent_variable <- "Sample_Group"
  model_result <- mean_permutation(family_test, sig.formula, tempDataFrame, independent_variable)

  prepared_data <- data_preparation(family_test,transformation,tempDataFrame, independent_variable, g_start, dototal, covariates, depth_analysis)
  tempDataFrame <- prepared_data$tempDataFrame
  independent_variable1stLevel <- prepared_data$independent_variable1stLevel
  independent_variable2ndLevel <- prepared_data$independent_variable2ndLevel

  Breusch_Pagan_pvalue <- NA
  cols <- colnames(tempDataFrame)
  iters <- length(cols)
  g <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = g_start:iters)
  else
    progress_bar <- ""

  to_export <- c("cols", "family_test", "covariates", "independent_variable", "tempDataFrame",
    "independent_variable1stLevel", "independent_variable2ndLevel",
    "key", "transformation","exact_pvalue","iters",
    "data_preparation","apply_stat_model_sig.formula","quantreg_permutation_model",
    "apply_stat_model_sig_formula", "data_distribution_info", "glm_model", "test_model", "Breusch_Pagan_pvalue",
    "progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_values","ssEnv","g_start")

 log_event("DEBUG: ", Sys.time(),  "Starting foreach withh: ", iters, " items")

 log_event("DEBUG: ", Sys.time(), " I'll perform:",iters - length(covariates)," tests." )
  result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind, .export = to_export) %dorng%
  # for(g in g_start:iters)
  {
    burdenValue <- cols[g]
    # update progress every 10%
    # ten_perc <- round(iters/100)
    if(ssEnv$showprogress)
      progress_bar(sprintf("doing genomic area: %s", stringr::str_pad(burdenValue, 10, pad = " ")))

    #
    if(!is.null(tempDataFrame[,burdenValue]) & length(unique(tempDataFrame[,burdenValue]))>=2){

      sig.formula <- apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
      bartlett_pvalue <- data_distribution_info(family_test, tempDataFrame, burdenValue, independent_variable)

      if( family_test=="binomial" | family_test=="poisson" | family_test=="gaussian")
        model_result <- glm_model(family_test, tempDataFrame, sig.formula )

      if(family_test=="wilcoxon" | family_test=="t.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" | family_test=="jsd"
        | family_test=="chisq.test" | family_test=="fisher.test")
        model_result <- test_model(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable )

      if(grepl("mean-permutation",family_test))
        model_result <- mean_permutation(family_test, sig.formula, tempDataFrame, independent_variable)

      if(grepl("quantile-permutation",family_test))
        model_result <- quantile_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable)

      if(grepl("spearman-permutation",family_test))
        model_result <- spearman_permutation(family_test, sig.formula, tempDataFrame, independent_variable)

      # Determine the null device for the current platform
      null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
      # Redirect output to the null device
      sink(null_device)
      if (grepl("quantreg-permutation", family_test))
        model_result <- quantreg_permutation_model(family_test, sig.formula, tempDataFrame, independent_variable)
      sink()
      # sink(file.path(ssEnv$session_folder,"session_output.log"), split = TRUE, append = TRUE)

      pvalue <- model_result$pvalue
      statistic_parameter <- model_result$statistic_parameter
      pvalueadjusted <- model_result$pvalue
      aic_value <- model_result$aic_value
      residuals <- model_result$residuals
      shapiro_pvalue <- model_result$shapiro_pvalue
      ci.lower <- model_result$ci.lower
      ci.upper <- model_result$ci.upper
      r_model <- model_result$r_model
      std.error <- model_result$std.error
      n_permutations <- model_result$n_permutations
      n_permutations_test <- model_result$n_permutations_test

      local_result <- data.frame("INDIPENDENT.VARIABLE" = independent_variable)
      local_result$MARKER <- as.character(key$MARKER)
      local_result$FIGURE <-  as.character(key$FIGURE)
      local_result$AREA <-  as.character(key$AREA)
      local_result$SUBAREA <-  as.character(key$SUBAREA)
      local_result$AREA_OF_TEST <- burdenValue
      local_result$PVALUE <- (model_result$pvalue)
      local_result$PVALUEADJ <- pvalueadjusted
      local_result$TEST <- "SINGLE_AREA"
      local_result$STATISTIC_PARAMETER <- if(exists("statistic_parameter")) statistic_parameter  else NA
      local_result$STD.ERROR <- if(exists("std.error")) std.error  else NA
      local_result$AIC <- if(exists("aic_value")) aic_value  else NA
      local_result$RESIDUALS.SUM <- if(exists("model_result")) (sum(model_result$residuals))  else NA
      local_result$FAMILY <- family_test
      local_result$R.PACK <- r_model
      local_result$transformation <- transformation
      local_result$COVARIATES <- paste0(covariatescollapse=" ")
      local_result$SHAPIRO.PVALUE <- shapiro_pvalue
      local_result$BREUSCH_PAGAN.PVALUE <- Breusch_Pagan_pvalue
      local_result$BARTLETT.PVALUE <- if(exists("bartlett_pvalue")) (bartlett_pvalue)  else NA
      local_result$RHO <- if(r_model=="stats_cor.test") model_result$statistic_parameter else NA
      local_result$CI.LOWER <- if(exists("ci.lower")) ci.lower else NA
      local_result$CI.UPPER <- if(exists("ci.upper")) ci.upper else NA
      local_result$N.PERMUTATIONS <- if(exists("n_permutations")) n_permutations else NA

      if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" &
         family_test!="kendall" & !grepl("quantreg-permutation", family_test)
         & family_test!="poisson"  & !grepl("mean-permutation", family_test))
      {
        independent_variableData1stLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable1stLevel,burdenValue])
        independent_variableData2ndLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable2ndLevel,burdenValue])

        local_result$CASE.LABEL <- if(exists("independent_variable1stLevel")) as.character(independent_variable1stLevel) else NA
        local_result$COUNT.CASE <- length(independent_variableData1stLevel)
        local_result$MEAN.CASE <- (mean(independent_variableData1stLevel))
        local_result$SD.CASE <- (stats::sd(independent_variableData1stLevel))
        local_result$CONTROL.LABEL <- as.character(independent_variable2ndLevel)
        local_result$COUNT.CONTROL <- length(stats::na.omit(independent_variableData2ndLevel))
        local_result$MEAN.CONTROL <- (mean(independent_variableData2ndLevel))
        local_result$SD.CONTROL <- (stats::sd(independent_variableData2ndLevel))
      } else
      {
        dependentVariableData <- as.numeric(stats::na.omit(tempDataFrame[!is.na(tempDataFrame[,independent_variable]),burdenValue]))
        independent_variableData <- as.numeric(stats::na.omit(tempDataFrame[  ,independent_variable]))

        if(sum(is.na(dependentVariableData)>0) | sum(is.na(independent_variableData)))
        {
          log_event("ERROR: ", Sys.time(), "The submitted data are not factorial or numeric.")
          stop()
        }
        local_result$CASE.LABEL <- NA
        local_result$COUNT.CASE <- NA
        local_result$MEAN.CASE <- NA
        local_result$SD.CASE <- NA
        local_result$CONTROL.LABEL <- NA
        local_result$COUNT.CONTROL <- NA
        local_result$MEAN.CONTROL <- NA
        local_result$SD.CONTROL <- NA
      }
      local_result
    }
  }

  #
  log_event("INFO: ", Sys.time(), " I performed:",iters," tests." )

  if(exists("result_temp") )
  {

    # result_temp <- unique(result_temp)
    result_temp <- result_temp %>% dplyr::distinct()

    if (!is.null(dim(result_temp)))
    {
      result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"))
      result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"))
    }
    return(result_temp)
  }
  return(NULL)
}
