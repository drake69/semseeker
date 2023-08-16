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

  boot_success <- if(is.null(arguments[["boot_success"]])) 0 else arguments$boot_success
  tests_count <- if(is.null(arguments[["tests_count"]])) 1 else arguments$tests_count
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
    "key", "transformation","quantreg_summary","iters", "boot_success", "tests_count",
    "data_preparation","apply_stat_model_sig.formula","quantreg_model",
    "apply_stat_model_sig_formula", "data_distribution_info", "glm_model", "test_model", "Breusch_Pagan_pvalue",
    "progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","beta_values","iqrTimes","ssEnv")
  # message("Starting foreach withh: ", iters, " items")

  message("INFO: ", Sys.time(), " I'll perform:",iters - length(covariates)," tests." )
  result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind, .export = to_export) %dorng%
  # for(g in g_start:iters)
  {
    burdenValue <- cols[g]
    if(ssEnv$showprogress)
      progress_bar(sprintf("genomic area: %s", stringr::str_pad( burdenValue, 20, side=c('left'), pad=' ')))
    if(!is.null(tempDataFrame[,burdenValue]) & length(unique(tempDataFrame[,burdenValue]))>2){

      sig.formula <- apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
      bartlett_pvalue <- data_distribution_info(family_test, tempDataFrame, burdenValue, independent_variable)

      if( family_test=="binomial" | family_test=="poisson" | family_test=="gaussian")
        model_result <- glm_model(family_test, tempDataFrame, sig.formula )

      if(family_test=="wilcoxon" | family_test=="t.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
        model_result <- test_model(family_test, tempDataFrame, sig.formula,burdenValue,independent_variable )

      # Determine the null device for the current platform
      null_device <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
      # Redirect output to the null device
      sink(null_device)

      # if (sink.number() != 0)
      #   sink(NULL)
      # sink(file.path(ssEnv$session_folder,"session_output.log"), split = FALSE, append = TRUE)
      # sink(tempfile())
      if (grepl("quantreg", family_test))
        model_result <- quantreg_model(family_test, sig.formula, tempDataFrame, independent_variable, boot_success, tests_count)
      sink()
      # sink(file.path(ssEnv$session_folder,"session_output.log"), split = TRUE, append = TRUE)

      pvalue <- model_result$pvalue
      beta_value <- model_result$beta_value
      pvalueadjusted <- model_result$pvalue
      aic_value <- model_result$aic_value
      residuals <- model_result$residuals
      shapiro_pvalue <- model_result$shapiro_pvalue
      ci.lower <- model_result$ci.lower
      ci.upper <- model_result$ci.upper
      ci.upper.adjusted <- model_result$ci.upper.adjusted
      ci.lower.adjusted <- model_result$ci.lower.adjusted
      r_model <- model_result$r_model
      std.error <- model_result$std.error
      n_permutations <- model_result$n_permutations

      if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" &
         family_test!="kendall" & !grepl("quantreg", family_test)
         & family_test!="poisson")
      {
        independent_variableData1stLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable1stLevel ,burdenValue])
        independent_variableData2ndLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable2ndLevel,burdenValue])
        local_result <- data.frame (
          "INDIPENDENT.VARIABLE"= independent_variable,
          "MARKER" = key$MARKER,
          "FIGURE" = key$FIGURE,
          "AREA" = key$AREA,
          "SUBAREA" = key$SUBAREA,
          "AREA_OF_TEST" = burdenValue,
          "PVALUE" = (model_result$pvalue),
          "PVALUEADJ" = pvalueadjusted,
          "TEST" = "SINGLE_AREA",
          "BETA" = if(exists("beta_value")) beta_value  else NA,
          "STD.ERROR" = if(exists("std.error")) std.error  else NA,
          "AIC" = if(exists("aic_value")) aic_value  else NA,
          "RESIDUALS.SUM" = if(exists("model_result")) (sum(model_result$residuals))  else NA,
          "FAMILY" = family_test,
          "R.PACK" = r_model,
          "transformation" = transformation,
          "COVARIATES" = paste0(covariates,collapse=" "),
          "SHAPIRO.PVALUE" = shapiro_pvalue,
          "BREUSCH-PAGAN.PVALUE" = Breusch_Pagan_pvalue,
          "BARTLETT.PVALUE" = if(exists("bartlett_pvalue")) (bartlett_pvalue)  else NA,
          "CASE.LABEL"= if(exists("independent_variable1stLevel")) as.character(independent_variable1stLevel) else NA,
          "COUNT.CASE"=length(independent_variableData1stLevel),
          "MEAN.CASE" = (mean(independent_variableData1stLevel)),
          "SD.CASE"=(stats::sd(independent_variableData1stLevel)),
          "CONTROL.LABEL" = as.character(independent_variable2ndLevel),
          "COUNT.CONTROL"=length(stats::na.omit(independent_variableData2ndLevel)),
          "MEAN.CONTROL"=(mean(independent_variableData2ndLevel)),
          "SD.CONTROL"= (stats::sd(independent_variableData2ndLevel)),
          "RHO"= if(r_model=="stats_cor.test") model_result$beta_value else NA,
          "CI.LOWER"= if(exists("ci.lower")) ci.lower else NA,
          "CI.UPPER"= if(exists("ci.upper")) ci.upper else NA,
          "CI.LOWER.ADJUSTED"=  if(exists("ci.lower.adjusted")) ci.lower.adjusted else NA,
          "CI.UPPER.ADJUSTED"=  if(exists("ci.upper.adjusted")) ci.upper.adjusted else NA,
          "N.PERMUTATIONS" =  NA
        )
      } else
      {
        dependentVariableData <- as.numeric(stats::na.omit(tempDataFrame[!is.na(tempDataFrame[,independent_variable]),burdenValue]))
        independent_variableData <- as.numeric(stats::na.omit(tempDataFrame[  ,independent_variable]))

        if(sum(is.na(dependentVariableData)>0) | sum(is.na(independent_variableData)))
        {
          message("ERROR: ", Sys.time(), "The submitted data are not factorial or numeric.")
          stop()
        }

        local_result <- data.frame (
          "INDIPENDENT.VARIABLE"= independent_variable,
          "MARKER" = key$MARKER,
          "FIGURE" = key$FIGURE,
          "AREA" = key$AREA,
          "SUBAREA" = key$SUBAREA,
          "AREA_OF_TEST" = burdenValue,
          "PVALUE" = (model_result$pvalue),
          "PVALUEADJ" = pvalueadjusted,
          "TEST" = "SINGLE_AREA",
          "BETA" = if(exists("beta_value")) beta_value  else NA,
          "STD.ERROR" = if(exists("std.error")) std.error  else NA,
          "AIC" = if(exists("aic_value")) aic_value  else NA,
          "RESIDUALS.SUM" = if(exists("model_result")) (sum(model_result$residuals))  else NA,
          "FAMILY" = family_test,
          "R.PACK" = r_model,
          "transformation" = transformation,
          "COVARIATES" = paste0(covariates,collapse=" "),
          "SHAPIRO.PVALUE" = shapiro_pvalue,
          "BREUSCH-PAGAN.PVALUE" = Breusch_Pagan_pvalue,
          "BARTLETT.PVALUE" = NA,
          "CASE.LABEL"= as.character(independent_variable),
          "COUNT.CASE"=length(independent_variableData),
          "MEAN.CASE" = NA,
          "SD.CASE"= NA,
          "CONTROL.LABEL" = "BURDEN.VALUE",
          "COUNT.CONTROL"=length(dependentVariableData),
          "MEAN.CONTROL"=(mean(dependentVariableData)),
          "SD.CONTROL"= (stats::sd(dependentVariableData)),
          "RHO"= if(r_model=="stats_cor.test")  (model_result$beta_value) else NA,
          "CI.LOWER"= if(exists("ci.lower")) ci.lower else NA,
          "CI.UPPER"= if(exists("ci.upper")) ci.upper else NA,
          "CI.LOWER.ADJUSTED"=  if(exists("ci.lower.adjusted")) ci.lower.adjusted else NA,
          "CI.UPPER.ADJUSTED"=  if(exists("ci.upper.adjusted")) ci.upper.adjusted else NA,
          "N.PERMUTATIONS" = if(exists("n_permutations")) n_permutations else NA
        )
      }
      # if(exists("result_temp"))
      #   result_temp <- rbind(result_temp, local_result)
      # else
      #   result_temp <- local_result

      local_result
    }
    else
    {
      # local_result <- data.frame (
      #   "INDIPENDENT.VARIABLE"= independent_variable,
      #   "MARKER" = key$MARKER,
      #   "FIGURE" = key$FIGURE,
      #   "AREA" = key$AREA,
      #   "SUBAREA" = key$SUBAREA,
      #   "AREA_OF_TEST" = burdenValue,
      #   "PVALUE" = NA,
      #   "PVALUEADJ" = NA,
      #   "TEST" = "SINGLE_AREA",
      #   "BETA" = NA,
      #   "STD.ERROR" = NA,
      #   "AIC" = NA,
      #   "RESIDUALS.SUM" = NA,
      # "FAMILY" = family_test,
      # "R.PACK" = r_model,
      #   "transformation" = transformation,
      #   "COVARIATES" = paste0(covariates,collapse=" "),
      #   "SHAPIRO.PVALUE" = NA,
      #   "BREUSCH-PAGAN.PVALUE" = NA,
      #   "BARTLETT.PVALUE" = NA,
      #   "CASE.LABEL"= NA,
      #   "COUNT.CASE"=NA,
      #   "MEAN.CASE" = NA,
      #   "SD.CASE"= NA,
      #   "CONTROL.LABEL" = "BURDEN.VALUE",
      #   "COUNT.CONTROL"= NA,
      #   "MEAN.CONTROL"= NA,
      #   "SD.CONTROL"= NA,
      #   "RHO"=  NA,
      #   "CI.LOWER"=  NA,
      #   "CI.UPPER"=  NA,
      #   "CI.LOWER.ADJUSTED"=  NA,
      #   "CI.UPPER.ADJUSTED"=  NA,
      #   "N.PERMUTATIONS" =  NA
      # )
    }
  }

  message("INFO: ", Sys.time(), " I performed:",iters - 2," tests." )

  if(exists("result_temp"))
  {
    result_temp <- unique(result_temp)
    result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"))
    result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"))
    return(result_temp)
  }
  return(NULL)
  if(ssEnv$showprogress)
    remove(progress_bar)
}
