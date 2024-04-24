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

  prepared_data <- data_preparation(family_test,transformation,tempDataFrame, independent_variable, g_start, dototal, covariates, depth_analysis)
  tempDataFrame <- prepared_data$tempDataFrame
  independent_variable1stLevel <- prepared_data$independent_variable1stLevel
  independent_variable2ndLevel <- prepared_data$independent_variable2ndLevel

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

  result_columns <- c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST", "CI.LOWER", "CI.UPPER", "PVALUE", "STATISTIC_PARAMETER", "AIC_VALUE", "RESIDUALS", "SHAPIRO_PVALUE", "R_MODEL", "STD.ERROR", "N_PERMUTATIONS", "N_PERMUTATIONS_TEST")
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  "Starting foreach withh: ", iters, " items")

  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " I'll perform:",iters - length(covariates)," tests." )

  result_temp <- foreach::foreach(g = g_start:iters, .combine =  plyr::rbind.fill, .export = to_export) %dorng%
  # for(g in g_start:iters)
  {
    burdenValue <- cols[g]
    if(ssEnv$showprogress)
      progress_bar(sprintf("doing genomic area: %s", stringr::str_pad(burdenValue, 10, pad = " ")))

    if(!is.null(tempDataFrame[,burdenValue]) & length(unique(tempDataFrame[,burdenValue]))>=2){

      sig.formula <- apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
      bartlett_pvalue <- data_distribution_info(family_test, tempDataFrame, burdenValue, independent_variable)
      model_result <- execute_model(family_test, tempDataFrame, sig.formula, burdenValue, independent_variable)

      local_result <- data.frame("INDIPENDENT.VARIABLE" = independent_variable)
      local_result$MARKER <- as.character(key$MARKER)
      local_result$FIGURE <-  as.character(key$FIGURE)
      local_result$AREA <-  as.character(key$AREA)
      local_result$SUBAREA <-  as.character(key$SUBAREA)
      local_result$AREA_OF_TEST <- burdenValue
      local_result$FAMILY <- family_test
      local_result$TEST <- "SINGLE_AREA"
      local_result$transformation <- transformation
      local_result$COVARIATES <- paste0(covariatescollapse=" ")

      if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" &
         family_test!="kendall" & !grepl("quantreg-permutation", family_test)
         & family_test!="poisson"  & !grepl("mean-permutation", family_test) &
          !grepl("polynomial", family_test) & !grepl("exp", family_test))
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
          log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "The submitted data are not factorial or numeric.")
          stop()
        }
      }

      colnames(model_result) <- toupper(colnames(model_result))
      colnames(local_result) <- toupper(colnames(local_result))
      local_result <- cbind(local_result, model_result)
      gc()
      local_result
    }
  }


  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I performed:",iters," tests." )

  gc()
  # & !is.null(result_temp)
  if(exists("result_temp") & !is.null(result_temp))
  {

    # result_temp <- unique(result_temp)
    result_temp <- result_temp %>% dplyr::distinct()

    if (!is.null(dim(result_temp)) )
    {
      if ("PVALUE" %in% colnames(result_temp))
      {
        result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"))
        result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"))
      }
    }
    return(result_temp)
  }
  return(NULL)
}
