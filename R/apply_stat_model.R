#' Title
#'
#' @param tempDataFrame data frame to apply association
#' @param g_start index of starting data
#' @param family_test family of test to run
#' @param covariates vector of covariates
#' @param key key to identify file to elaborate
#' @param transformation transformation to apply to covariates, burden and independent variable
#' @param dototal do a total per area
#' @param logFolder where to save log file
#' @param independent_variable independent variable name
#' @param depth_analysis depth's analysis
#' @param envir object environment
#' @param ... extra parameters
#'
#' @importFrom doRNG %dorng%
#'
apply_stat_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation, dototal, logFolder,
                             independent_variable, depth_analysis=3, envir , ...)
{
  arguments <- list(...)
  boot_success <- if(is.null(arguments[["boot_success"]])) 0 else arguments$boot_success
  tests_count <- if(is.null(arguments[["tests_count"]])) 1 else arguments$tests_count

  transformation <- as.character(transformation)
  originalDataFrame <- tempDataFrame

  if(is.factor(tempDataFrame[, independent_variable]))
  {
    independent_variable1stLevel <- levels(tempDataFrame[, independent_variable])[1]
    independent_variable2ndLevel <- levels(tempDataFrame[, independent_variable])[2]
  }

  df_head <- tempDataFrame[,1:(g_start-1)]
  burden_values <- sapply(tempDataFrame[,g_start:ncol(tempDataFrame)], as.numeric)


  df_colnames <- colnames(tempDataFrame)
  if( !is.null(dim(burden_values))  & dototal) {
    sum_area <- apply(burden_values, 1, sum)
    if(depth_analysis==2)
    {
      #select just column of independent variables, remove columns burden value, preserve only total
      df_colnames <- c(df_colnames[!(df_colnames %in% colnames(burden_values))],"TOTAL")
      burden_values <- data.frame("TOTAL"=sum_area)
    }
    else
    {
      burden_values <- data.frame(burden_values,"TOTAL"=sum_area)
      df_colnames <- c(df_colnames,"TOTAL")
    }
  }

  if(family_test != "poisson")
    burden_values <- burden_values + 0.001

  transformation <- as.character(transformation)
  if(is.null(transformation) | length(transformation)==0 | is.na(transformation))
    transformation <- "none"


  burden_values <- as.data.frame(burden_values)
  df_values_orig <- burden_values
  try(
    {
      burden_values = switch(
        as.character(transformation),
        "scale" = if(ncol(burden_values)>1) as.data.frame(apply(burden_values,2,scale)) else scale(burden_values),
        "log" = log(burden_values),
        "log2" = log2(burden_values),
        "log10"= log10(burden_values),
        "exp" = exp(burden_values),
        "none" = burden_values,
        burden_values
      )
    }
  )
  if(grepl("quantile", transformation))
  {
    qq <- as.numeric(unlist(strsplit(transformation,"\\_"))[2])
    burden_values <- as.data.frame(apply(burden_values,2,function(x){
      if(length(unique(x))>=qq)
        as.numeric(dplyr::ntile(x, n=qq))
      else
        rep(0,length(x))
      }))
  }
  burden_values <- as.data.frame(burden_values)

  if(setequal(burden_values,df_values_orig) & transformation !="none")
    transformation <- paste0("NA_", transformation, sep="")


  if(family_test!="binomial" & family_test!="wilcoxon" & family_test!="t.test" & family_test!="poisson")
  {
    variable_to_transform <- independent_variable
    if(length(covariates)>0)
    {
      variable_to_transform <- c(independent_variable,covariates)
      independent_variableValues <-as.data.frame(apply(tempDataFrame[,variable_to_transform] ,2, as.numeric))
    }
    else
    {
      independent_variableValues <-as.data.frame(as.numeric(tempDataFrame[,variable_to_transform]))
    }
    independent_variableValuesOrig <- independent_variableValues
    try(
      {
        independent_variableValues = switch(
          as.character(transformation),
          "scale" = if(ncol(independent_variableValues)>1) as.data.frame(apply(independent_variableValues,2,scale)) else scale(independent_variableValues),
          "log" = log(independent_variableValues),
          "log2" = log2(independent_variableValues),
          "log10"= log10(independent_variableValues),
          "exp" = exp(independent_variableValues),
          "none" = independent_variableValues,
          independent_variableValues
        )
      }
    )

    # if(grepl("quantile", transformation))
    # {
    #   qq <- unlist(strsplit(transformation,"_")[2])
    #   df_values_temp <- as.data.frame(apply( burden_values,2,function(x) dplyr::ntile(x, n=qq)))
    #   colnames(df_values_temp) <- colnames(burden_values)
    # }

    if(setequal(burden_values,df_values_orig) & transformation !="none")
      transformation <- paste0("NA_", transformation, sep="")
    else
      tempDataFrame[, variable_to_transform] <- independent_variableValues
  }


  tempDataFrame <- data.frame(df_head, burden_values)
  if(ncol(tempDataFrame)!=length(df_colnames))
    browser()
  colnames(tempDataFrame) <- df_colnames
  cols <- colnames(tempDataFrame)
  iters <- length(cols)


  # after the transformation some data could be missed
  lost_cols <- colSums(apply(tempDataFrame,2,is.nan))!=0
  lostDataFrame <-  colnames(tempDataFrame)[lost_cols]
  if(sum(lost_cols)!=0)
    utils::write.csv2(lostDataFrame, file.path(envir$logFolder,paste("lost_data_",transformation,"_",stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),".log", sep="")))

  #  we want to preserve the NA in the indipendent variables to be removed by the models
  tempDataFrame[apply(tempDataFrame,2,is.nan)] <- 0

  g <- 0
  to_export <- c("cols", "family_test", "covariates", "independent_variable", "tempDataFrame",
                 "independent_variable1stLevel", "independent_variable2ndLevel",
                 "key", "transformation","quantreg_summary","iters", "boot_success", "tests_count")

  # message("Starting foreach withh: ", iters, " items")
  message(Sys.time()," I'll perform:",iters," tests." )
  result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind, .export = to_export) %dorng%
  # for(g in g_start:iters)
  {
    #g <- 2
    burdenValue <- cols[g]
    # message(g)
    if(!is.null(tempDataFrame[,burdenValue]) & length(unique(tempDataFrame[,burdenValue]))>2){

      if(family_test=="wilcoxon" | family_test=="t.test")
      {
        covariates_model <- independent_variable
        sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
      }

      if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
      {
        covariates_model <- independent_variable
        sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
      }

      if (family_test=="binomial")
      {
        # inversion of roles for variable
        if(is.null(covariates) || length(covariates)==0)
        {
          covariates_model <- burdenValue
        } else
        {
          covariates_model <- paste0(paste0(c(burdenValue, covariates),collapse="+", sep=""))
        }
        sig.formula <- stats::as.formula(paste0(independent_variable,"~", covariates_model, sep=""))
      }

      if(family_test=="gaussian" | family_test=="poisson" | grepl("quantreg", family_test))
      {
        if(is.null(covariates) || length(covariates)==0)
          covariates_model <- independent_variable
        else
          covariates_model <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
        sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
      }

      if(family_test=="binomial" | family_test=="wilcoxon" | family_test=="t.test")
      {
        bartlett_pvalue <- stats::bartlett.test( stats::as.formula(paste0(burdenValue,"~", independent_variable, sep="")),
                                                 data= as.data.frame(tempDataFrame) )
      }

      if(family_test=="gaussian" | family_test=="spearman" | family_test=="kendall" | family_test=="pearson")
      {
        localDataFrame <- data.frame("depVar"=tempDataFrame[, burdenValue],"indepVar"=1 )
        localDataFrame <- rbind( localDataFrame,  data.frame("depVar"=tempDataFrame[, independent_variable],"indepVar"=2 ))
        bartlett_pvalue <- stats::bartlett.test( stats::as.formula("depVar ~ indepVar"), data= localDataFrame )
      }


      Breusch_Pagan_pvalue <- NA
      if( family_test=="binomial" | family_test=="poisson")
      {
        result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
        pvalue <- summary(result_glm )$coeff[-1, 4][1]
        beta_value <- (summary(result_glm )$coeff[-1, 1][1])
        aic_value <- (result_glm$aic)
        residuals <-  result_glm$resid
        #calculate shapiro of working residuals
        shapiro_pvalue <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
        # Breusch_Pagan_pvalue <- lmtest::bptest( data=residuals )$p.value
      } else if(family_test=="gaussian")
      {

        # sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
        # model <- stats::lm(formula = sig.formula, data = as.data.frame(tempDataFrame))
        # residuals <-  model$residuals
        # shapiro_pvalue <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
        # Breusch_Pagan_pvalue <- lmtest::bptest(model)$p.value

        result_glm  <- stats::lm( formula = sig.formula, data = as.data.frame(tempDataFrame))
        pvalue <- summary(result_glm )$coeff[-1, 4][1]
        beta_value <- (summary(result_glm )$coeff[-1, 1][1])
        aic_value <- 0 #(result_glm$aic)
        residuals <-  result_glm$residuals
        #calculate shapiro of working residuals
        shapiro_pvalue <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
        # Breusch_Pagan_pvalue <- lmtest::bptest(data=residuals )$p.value
      }
      else
        #calculate shapiro of burden values
        shapiro_pvalue <- if(length(tempDataFrame[,burdenValue])>3 & length(unique(tempDataFrame[,burdenValue]))>3) (stats::shapiro.test(tempDataFrame[,burdenValue])$p.value) else NA

      if(grepl("quantreg", family_test))
      {
        lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
        quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
        if(length(quantreg_params)<5)
        {
          if(length(quantreg_params)<3)
          {
            message("Nothing to do! Not enough parameter for quantile regression!.")
            return(NULL)
          }
          else
          {
            # Define function to compute p-value and beta regression coefficient
            compute_beta <- function(sig.formula, tau, dataFrame) {
              fit <- quantreg::rq(formula =  sig.formula,data = as.data.frame(dataFrame),  tau = tau)
              coef <-as.data.frame(summary(fit, se = "boot")$coefficients)[2,"Value"]
              pval <- as.data.frame(summary(fit, se = "boot")$coefficients)[2, "Pr(>|t|)"]
              return(list(beta = coef, pval = pval))
            }

            tau = as.numeric(quantreg_params[2])
            n_permutations <- as.numeric(quantreg_params[3])
            # Compute beta and p-value for n_permutations replications
            results <- replicate(n_permutations, compute_beta(sig.formula, tau, dataFrame=tempDataFrame))
            # Compute average beta and p-value
            beta_value <- mean(unlist(t(results)[,"beta"]))
            pvalue <- mean(unlist(t(results)[,"pval"]))

            # tau = as.numeric(quantreg_params[2])
            # result_quantreg <- quantreg::rq(formula = sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, method="sfn", model = T)
            # pvalue <- summary(result_quantreg )[["tTable"]][2,5]
            # beta_value <- (summary(result_quantreg)$coeff[-1, 1][1])
            # aic_value <- (result_quantreg$aic)
            # residuals <-  result_quantreg$resid
            #calculate shapiro of working residuals
            # shapiro_pvalue <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
          }
        }
        else
        {
          n_permutations_test <- as.numeric(quantreg_params[3])
          n_permutations <- as.numeric(quantreg_params[4])
          tau <- as.numeric(quantreg_params[2])
          conf.level <- as.numeric(quantreg_params[5])
          model.x <-  suppressMessages(lqmm::lqm( formula = sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control))
          if(n_permutations > n_permutations_test)
          {
            model.x.boot <- suppressMessages(lqmm::boot(model.x, R = n_permutations_test))
            beta_value <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
            tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
            colnames(tt) <- colnames(model.x.boot)
            boot_vector <- stats::na.omit(tt[,independent_variable])
            boot.bca <- quantreg_summary(boot_vector, beta_value, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable, lqm_control = lqm_control, conf.level = conf.level)
          }
          ci.lower.adjusted <- NA
          ci.upper.adjusted <- NA

          if(boot.bca[1] <0 & boot.bca[2]>0 )
          {
            n_permutations <- n_permutations_test
          }
          else
          {
            model.x.boot <- suppressMessages(lqmm::boot(model.x, R = n_permutations))
            beta_value <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
            tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
            colnames(tt) <- colnames(model.x.boot)
            boot_vector <- stats::na.omit(tt[,independent_variable])
            boot.bca <- quantreg_summary(boot_vector, beta_value, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable, lqm_control = lqm_control, conf.level = conf.level)
            boot.bca.adjusted <- quantreg_summary(boot_vector, beta_value, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable, lqm_control = lqm_control,  boot_success = boot_success, tests_count=tests_count, conf.level = conf.level)
            ci.lower.adjusted <-  boot.bca.adjusted[1]
            ci.upper.adjusted <- boot.bca.adjusted[2]
          }

          ci.lower <- boot.bca[1]
          ci.upper <- boot.bca[2]
          pvalue <- boot.bca[3]
        }
      }

      if(family_test=="wilcoxon")
      {
        result_w  <- suppressWarnings(stats::wilcox.test(formula= sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE))
        pvalue <- result_w$p.value
      }

      if(family_test=="t.test")
      {
        result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame))
        pvalue <- result_w$p.value
      }

      if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
      {

        result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
        pvalue <- result_cor$p.value
      }


      pvalueadjusted <- pvalue
      # beta <- exp(summary( result_glm )$coeff[-1, 1][1])
      # result_temp <-

      if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" &
         family_test!="kendall" & !grepl("quantreg", family_test)
         & family_test!="poisson")
      {
        independent_variableData1stLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable1stLevel ,burdenValue])
        independent_variableData2ndLevel <- stats::na.omit(tempDataFrame[tempDataFrame[, independent_variable]==independent_variable2ndLevel,burdenValue])
        local_result <- data.frame (
          "INDIPENDENT.VARIABLE"= independent_variable,
          "ANOMALY" = key$ANOMALY,
          "FIGURE" = key$FIGURE,
          "GROUP" = key$GROUP,
          "SUBGROUP" = key$SUBGROUP,
          "AREA_OF_TEST" = burdenValue,
          "PVALUE" = (pvalue),
          "PVALUEADJ" = pvalueadjusted,
          "TEST" = "SINGLE_AREA",
          "BETA" = if(exists("beta_value")) beta_value  else NA,
          "AIC" = if(exists("aic_value")) aic_value  else NA,
          "RESIDUALS.SUM" = if(exists("result_glm")) (sum(result_glm$residuals))  else NA,
          "FAMILY" = family_test,
          "transformation" = transformation,
          "COVARIATES" = paste0(covariates,collapse=" "),
          "SHAPIRO.PVALUE" = shapiro_pvalue,
          "BREUSCH-PAGAN.PVALUE" = Breusch_Pagan_pvalue,
          "BARTLETT.PVALUE" = if(exists("bartlett_pvalue")) (bartlett_pvalue$p.value)  else NA,
          "CASE.LABEL"= if(exists("independent_variable1stLevel")) as.character(independent_variable1stLevel) else NA,
          "COUNT.CASE"=length(independent_variableData1stLevel),
          "MEAN.CASE" = (mean(independent_variableData1stLevel)),
          "SD.CASE"=(stats::sd(independent_variableData1stLevel)),
          "CONTROL.LABEL" = as.character(independent_variable2ndLevel),
          "COUNT.CONTROL"=length(stats::na.omit(independent_variableData2ndLevel)),
          "MEAN.CONTROL"=(mean(independent_variableData2ndLevel)),
          "SD.CONTROL"= (stats::sd(independent_variableData2ndLevel)),
          "RHO"= if(exists("result_cor"))  (result_cor$estimate) else NA,
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
        local_result <- data.frame (
          "INDIPENDENT.VARIABLE"= independent_variable,
          "ANOMALY" = key$ANOMALY,
          "FIGURE" = key$FIGURE,
          "GROUP" = key$GROUP,
          "SUBGROUP" = key$SUBGROUP,
          "AREA_OF_TEST" = burdenValue,
          "PVALUE" = (pvalue),
          "PVALUEADJ" = pvalueadjusted,
          "TEST" = "SINGLE_AREA",
          "BETA" = if(exists("beta_value")) beta_value  else NA,
          "AIC" = if(exists("aic_value")) aic_value  else NA,
          "RESIDUALS.SUM" = if(exists("result_glm")) (sum(result_glm$residuals))  else NA,
          "FAMILY" = family_test,
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
          "RHO"= if(exists("result_cor"))  (result_cor$estimate) else NA,
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
      local_result <- data.frame (
        "INDIPENDENT.VARIABLE"= independent_variable,
        "ANOMALY" = key$ANOMALY,
        "FIGURE" = key$FIGURE,
        "GROUP" = key$GROUP,
        "SUBGROUP" = key$SUBGROUP,
        "AREA_OF_TEST" = burdenValue,
        "PVALUE" = NA,
        "PVALUEADJ" = NA,
        "TEST" = "SINGLE_AREA",
        "BETA" = NA,
        "AIC" = NA,
        "RESIDUALS.SUM" = NA,
        "FAMILY" = family_test,
        "transformation" = transformation,
        "COVARIATES" = paste0(covariates,collapse=" "),
        "SHAPIRO.PVALUE" = NA,
        "BREUSCH-PAGAN.PVALUE" = NA,
        "BARTLETT.PVALUE" = NA,
        "CASE.LABEL"= NA,
        "COUNT.CASE"=NA,
        "MEAN.CASE" = NA,
        "SD.CASE"= NA,
        "CONTROL.LABEL" = "BURDEN.VALUE",
        "COUNT.CONTROL"= NA,
        "MEAN.CONTROL"= NA,
        "SD.CONTROL"= NA,
        "RHO"=  NA,
        "CI.LOWER"=  NA,
        "CI.UPPER"=  NA,
        "CI.LOWER.ADJUSTED"=  NA,
        "CI.UPPER.ADJUSTED"=  NA,
        "N.PERMUTATIONS" =  NA
      )
    }
  }

  message(Sys.time()," I performed:",iters," tests." )

  if(exists("result_temp"))
  {
    result_temp <- unique(result_temp)
    result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"))
    result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- (stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"))
    return(result_temp)
  }
  return(NULL)
}
