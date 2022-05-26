apply_stat_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation, dototal, logFolder, independent_variable, depthAnalysis=3)
{
  # #browser()
  # parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c())
  transformation <- as.character(transformation)
  #browser()

  if(is.factor(tempDataFrame[, independent_variable]))
  {
    independent_variable1stLevel <- levels(tempDataFrame[, independent_variable])[1]
    independent_variable2ndLevel <- levels(tempDataFrame[, independent_variable])[2]
  }

  # if(independent_variable=="Sample_Group")
  # {
  #   tempDataFrame[, independent_variable][tempDataFrame[, independent_variable]=="Control"] <- FALSE
  #   tempDataFrame[, independent_variable][tempDataFrame[, independent_variable]=="Case"] <- TRUE
  #   # tempDataFrame[, independent_variable] <- as.numeric(tempDataFrame[, independent_variable])
  # }


  df_head <- tempDataFrame[,1:(g_start-1)]
  df_values <- sapply(tempDataFrame[,g_start:ncol(tempDataFrame)], as.numeric)

  if( ncol(tempDataFrame) - g_start  > 1) {
    df_values <- df_values[, colSums(df_values) > 0]
  }

  df_colnames <- colnames(tempDataFrame)
  if( !is.null(dim(df_values))  & dototal) {
    sum_area <- apply(df_values, 1, sum)
    df_values <- data.frame(df_values,"TOTAL"=sum_area)
    df_colnames <- c(df_colnames,"TOTAL")
    if(depthAnalysis==2)
    {
      # #browser()
      df_values <- as.numeric(f_values[, colnames(df_values) %in% c( independent_variable, "TOTAL") ])
      df_colnames <- c( independent_variable, "TOTAL")
    }
  }

  # #browser()
  if(family_test != "poisson")
    df_values <- df_values + 0.001

  transformation <- as.character(transformation)
  if(is.null(transformation) | length(transformation)==0 | is.na(transformation))
    transformation <- "none"

  df_values_orig <- df_values
  try(
    {
      df_values = switch(
        as.character(transformation),
        "scale" = scale(df_values),
        "log" = log(df_values),
        "log2" = log2(df_values),
        "log10"= log10(df_values),
        "exp" = exp(df_values),
        "johnson" =  Johnson::RE.Johnson(df_values)$transformed,
        "none" = df_values,
        df_values
      )
    }
  )

  if(grepl("quantile", transformation))
  {
    qq <- strsplit(transformation,"_")[2]
    df_values_temp <- as.data.frame(apply( df_values,2,function(x) ggplot2::cut_number(x, n=qq)))
    colnames(df_values_temp) <- colnames(df_values)
  }

  if(setequal(df_values,df_values_orig) & transformation !="none")
    transformation <- paste0("NA_", transformation, sep="")


  if(family_test!="binomial" & family_test!="wilcoxon" & family_test!="t.test" & family_test!="poisson")
  {
    independent_variableValues <- as.numeric(tempDataFrame[, independent_variable])
    independent_variableValuesOrig <- as.numeric(tempDataFrame[, independent_variable])
    try(
      {
        independent_variableValues = switch(
          as.character(transformation),
          "scale" = scale(independent_variableValues),
          "log" = log(independent_variableValues),
          "log2" = log2(independent_variableValues),
          "log10"= log10(independent_variableValues),
          "exp" = exp(independent_variableValues),
          "johnson" =  Johnson::RE.Johnson(independent_variableValues)$transformed,
          "none" = independent_variableValues,
          independent_variableValues
        )
      }
    )

    if(setequal(df_values,df_values_orig) & transformation !="none")
      transformation <- paste0("NA_", transformation, sep="")
    else
    {
      tempDataFrame[, independent_variable] <- independent_variableValues
      }
  }

  tempDataFrame <- data.frame(df_head, df_values)
  # #browser()
  colnames(tempDataFrame) <- df_colnames
  cols <- colnames(tempDataFrame)
  iters <- length(cols)

  g <- 0
  to_export <- c("cols", "family_test", "covariates", "independent_variable", "tempDataFrame", "independent_variable1stLevel", "independent_variable2ndLevel",
                 "key", "transformation")
  result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind, .export = to_export) %dorng%
  # for (g in g_start:iters)
  {
    #g <- 2
    burdenValue <- cols[g]
    if (family_test=="poisson")
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- independent_variable
      } else
      {
        covariates_model <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
      }
      sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
    }

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

    if (family_test=="gaussian")
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- independent_variable
      } else
      {
        covariates_model <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
      }
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

    # #browser()
    if(family_test=="binomial" | family_test=="poisson" | family_test=="wilcoxon" | family_test=="t.test")
    {
      bartlett_pvalue <- stats::bartlett.test( stats::as.formula(paste0(burdenValue,"~", independent_variable, sep="")), data= as.data.frame(tempDataFrame) )
    }

    if(family_test=="gaussian" | family_test=="spearman" | family_test=="kendall" | family_test=="pearson")
    {
      localDataFrame <- data.frame("depVar"=tempDataFrame[, burdenValue],"indepVar"=1 )
      localDataFrame <- rbind( localDataFrame,  data.frame("depVar"=tempDataFrame[, independent_variable],"indepVar"=2 ))
      bartlett_pvalue <- stats::bartlett.test( stats::as.formula("depVar ~ indepVar"), data= localDataFrame )
    }

    shapiro_pvalue <- stats::shapiro.test(tempDataFrame[,burdenValue])

    if(family_test=="gaussian" | family_test=="binomial" | family_test=="poisson")
    {
      result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
      pvalue <- summary(result_glm )$coeff[-1, 4][1]
    }

    if(grepl("quantreg", family_test))
    {
      #browser()
      lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
      n_runs=strsplit(family_test,"_")[3]
      tau=strsplit(family_test,"_")[2]
      model.x <-  lqmm::lqm(sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control)
      model.x.boot <- lqmm::boot(model.x, R = n_runs)

      beta_full <- summary(model.x.boot)["independent_variable","Value"]
      tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
      colnames(tt) <- colnames(model.x.boot)

      # #####browser()
      boot_vector <- stats::na.omit(tt[,"independent_variable"])
      # boot.bca <- coxed::bca(boot_vector)
      boot.bca <- bca_pvalue_for_lqm(estimate = beta_full, boot_vector = boot_vector, model = model.x,working_data = tempDataFrame)
      boot_result <- data.frame("Value"=beta_full,"Bias"="","Std. Error"="","Lower bound"=boot.bca[1],"Upper bound"=boot.bca[2],"Pr(>|t|)"=as.numeric(boot.bca[3]))

      # result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
      pvalue <- boot_result[3]
    }

    if(family_test=="wilcoxon")
    {
      result_w  <- stats::wilcox.test(sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE)
      pvalue <- result_w$p.value
    }

    if(family_test=="t.test")
    {
      result_w  <-stats::t.test(sig.formula, data = as.data.frame(tempDataFrame))
      pvalue <- result_w$p.value
    }

    if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
    {

      result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
      pvalue <- result_cor$p.value
    }


    pvalueadjusted <- pvalue
    # #browser()
    # beta <- exp(summary( result_glm )$coeff[-1, 1][1])
    # result_temp <-

    if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" & family_test!="kendall")
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
        "PVALUE" = round(pvalue,3),
        "PVALUEADJ" = pvalueadjusted,
        "TEST" = "SINGLE_AREA",
        "BETA" = if(exists("result_glm")) round(summary(result_glm )$coeff[-1, 1][1],3)  else NA,
        "AIC" = if(exists("result_glm")) round(result_glm$aic,3)  else NA,
        "RESIDUALS.SUM" = if(exists("result_glm")) round(sum(result_glm$residuals),3)  else NA,
        "FAMILY" = family_test,
        "transformation" = transformation,
        "COVARIATES" = paste0(covariates,collapse=" "),
        "SHAPIRO.PVALUE" = round(shapiro_pvalue$p.value,3),
        "BARTLETT.PVALUE" = if(exists("bartlett_pvalue")) round(bartlett_pvalue$p.value,3)  else NA,
        "CASE.LABEL"= as.character(independent_variable1stLevel),
        "COUNT.CASE"=length(independent_variableData1stLevel),
        "MEAN.CASE" = round(mean(independent_variableData1stLevel),3),
        "SD.CASE"=round(stats::sd(independent_variableData1stLevel),3),
        "CONTROL.LABEL" = as.character(independent_variable2ndLevel),
        "COUNT.CONTROL"=length(stats::na.omit(independent_variableData2ndLevel)),
        "MEAN.CONTROL"=round(mean(independent_variableData2ndLevel),3),
        "SD.CONTROL"= round(stats::sd(independent_variableData2ndLevel),3),
        "RHO"= if(exists("result_cor"))  round(result_cor$estimate,4) else NA
      )
    } else
    {
      # #browser()
      dependentVariableData <- as.numeric(stats::na.omit(tempDataFrame[!is.na(tempDataFrame[,independent_variable]),burdenValue]))
      independent_variableData <- as.numeric(stats::na.omit(tempDataFrame[  ,independent_variable]))
      local_result <- data.frame (
        "INDIPENDENT.VARIABLE"= independent_variable,
        "ANOMALY" = key$ANOMALY,
        "FIGURE" = key$FIGURE,
        "GROUP" = key$GROUP,
        "SUBGROUP" = key$SUBGROUP,
        "AREA_OF_TEST" = burdenValue,
        "PVALUE" = round(pvalue,3),
        "PVALUEADJ" = pvalueadjusted,
        "TEST" = "SINGLE_AREA",
        "BETA" = if(exists("result_glm")) round(summary( result_glm )$coeff[-1, 1][1],3) else NA,
        "AIC" = if(exists("result_glm")) round(result_glm$aic,3)  else NA,
        "RESIDUALS.SUM" = if(exists("result_glm")) round(sum(result_glm$residuals),3)  else NA,
        "FAMILY" = family_test,
        "transformation" = transformation,
        "COVARIATES" = paste0(covariates,collapse=" "),
        "SHAPIRO.PVALUE" = round(shapiro_pvalue$p.value,3),
        "BARTLETT.PVALUE" = NA,
        "CASE.LABEL"= as.character(independent_variable),
        "COUNT.CASE"=length(independent_variableData),
        "MEAN.CASE" = round(mean(independent_variableData),3),
        "SD.CASE"= round(stats::sd(independent_variableData ),3),
        "CONTROL.LABEL" = "BURDEN.VALUE",
        "COUNT.CONTROL"=length(dependentVariableData),
        "MEAN.CONTROL"=round(mean(dependentVariableData),3),
        "SD.CONTROL"= round(stats::sd(dependentVariableData),3),
        "RHO"= if(exists("result_cor"))  round(result_cor$estimate,4) else NA
      )
    }
    # if (exists("result_temp"))
    # {
    #   result_temp <- rbind(result_temp, local_result)
    # }
    # else
    # {
    #   result_temp <- local_result
    # }
    local_result
  }

  result_temp <- unique(result_temp)


  result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- round(stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"),3)
  result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- round(stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"),3)

  gc()
  return(result_temp)
}
