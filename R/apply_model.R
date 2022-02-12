
apply_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation, dototal, logFolder, independentVariable, depthAnalysis=3)
{

  # browser()
  # parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c())

  if(independentVariable=="Sample_Group")
  {
    tempDataFrame[, independentVariable][tempDataFrame[, independentVariable]=="Control"] <- FALSE
    tempDataFrame[, independentVariable][tempDataFrame[, independentVariable]=="Case"] <- TRUE
    # tempDataFrame[, independentVariable] <- as.numeric(tempDataFrame[, independentVariable])
  }


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
      # browser()
      df_values <- df_values[, colnames(df_values) %in% c( independentVariable, "TOTAL") ]
      df_colnames <- c( independentVariable, "TOTAL")
    }
  }

  # browser()
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


  if(setequal(df_values,df_values_orig) & transformation !="none")
    transformation <- paste0("NA_", transformation, sep="")


  if(family_test!="binomial" & family_test!="wilcoxon" & family_test!="t.test" & family_test!="poisson")
  {
    independentVariableValues <- tempDataFrame[, independentVariable]
    independentVariableValuesOrig <- tempDataFrame[, independentVariable]
    try(
      {
        independentVariableValues = switch(
          as.character(transformation),
          "scale" = scale(independentVariableValues),
          "log" = log(independentVariableValues),
          "log2" = log2(independentVariableValues),
          "log10"= log10(independentVariableValues),
          "exp" = exp(independentVariableValues),
          "johnson" =  Johnson::RE.Johnson(independentVariableValues)$transformed,
          "none" = independentVariableValues,
          independentVariableValues
        )
      }
    )

    if(setequal(df_values,df_values_orig) & transformation !="none")
      transformation <- paste0("NA_", transformation, sep="")
    else
    {
      tempDataFrame[, independentVariable] <- independentVariableValues
      }
  }

  tempDataFrame <- data.frame(df_head, df_values)
  # browser()
  colnames(tempDataFrame) <- df_colnames
  cols <- colnames(tempDataFrame)
  iters <- length(cols)

  # result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind) %dopar%
  for (g in g_start:iters)
  {
    #g <- 2
    burdenValue <- cols[g]
    if (family_test=="poisson")
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- independentVariable
      } else
      {
        covariates_model <- paste0(paste0(c(independentVariable, covariates),collapse="+", sep=""))
      }
      sig.formula  <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
    }

    if(family_test=="wilcoxon" | family_test=="t.test")
    {
      covariates_model <- independentVariable
      sig.formul <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
    }

    if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
    {
      covariates_model <- independentVariable
      sig.formul <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
    }

    if (family_test=="gaussian")
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- independentVariable
      } else
      {
        covariates_model <- paste0(paste0(c(independentVariable, covariates),collapse="+", sep=""))
      }
      sig.formul <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
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
      sig.formul <- stats::as.formula(paste0(independentVariable,"~", covariates_model, sep=""))
    }

    # browser()
    if(family_test=="binomial" | family_test=="poisson" | family_test=="wilcoxon" | family_test=="t.test")
    {
      bartlett_pvalue <- stats::bartlett.test( stats::as.formula(paste0(burdenValue,"~", independentVariable, sep="")), data= as.data.frame(tempDataFrame) )
    }

    if(family_test=="gaussian" | family_test=="spearman" | family_test=="kendall" | family_test=="pearson")
    {
      localDataFrame <- data.frame("depVar"=tempDataFrame[, burdenValue],"indepVar"=1 )
      localDataFrame <- rbind( localDataFrame,  data.frame("depVar"=tempDataFrame[, independentVariable],"indepVar"=2 ))
      bartlett_pvalue <- stats::bartlett.test( stats::as.formula("depVar ~ indepVar"), data= localDataFrame )
    }

    shapiro_pvalue <- stats::shapiro.test(tempDataFrame[,burdenValue])

    if(family_test=="gaussian" | family_test=="binomial" | family_test=="poisson")
    {
      result.glm  <- stats::glm( sig.formula, family = family_test, data = as.data.frame(tempDataFrame))
      pvalue <- summary(result.glm )$coeff[-1, 4][1]
    }

    if(family_test=="wilcoxon")
    {
      result.w  <- stats::wilcox.test(sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE)
      pvalue <- result.w$p.value
    }

    if(family_test=="t.test")
    {
      result.w  <- stats::t.test(sig.formula, data = as.data.frame(tempDataFrame))
      pvalue <- result.w$p.value
    }

    if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
    {
      result.cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independentVariable]), method = family_test)
      pvalue <- result.cor$p.value
    }


    pvalueadjusted <- pvalue
    # browser()
    # beta <- exp(summary( result.glm )$coeff[-1, 1][1])
    # result_temp <-

    if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" & family_test!="kendall")
    {
      independentVariableData <- stats::na.omit(tempDataFrame[tempDataFrame[, independentVariable]==TRUE,burdenValue])
      dependentVariableData <- stats::na.omit(tempDataFrame[tempDataFrame[, independentVariable]==FALSE,burdenValue])
      local_result <- data.frame (
        "INDIPENDENT.VARIABLE"= independentVariable,
        "ANOMALY" = key$ANOMALY,
        "FIGURE" = key$FIGURE,
        "GROUP" = key$GROUP,
        "SUBGROUP" = key$SUBGROUP,
        "AREA_OF_TEST" = burdenValue,
        "PVALUE" = round(pvalue,3),
        "PVALUEADJ" = pvalueadjusted,
        "TEST" = "SINGLE_AREA",
        "BETA" = if(exists("result.glm")) round(summary(result.glm )$coeff[-1, 1][1],3)  else NA,
        "AIC" = if(exists("result.glm")) round(result.glm$aic,3)  else NA,
        "RESIDUALS.SUM" = if(exists("result.glm")) round(sum(result.glm$residuals),3)  else NA,
        "FAMILY" = family_test,
        "transformation" = transformation,
        "COVARIATES" = paste0(covariates,collapse=" "),
        "SHAPIRO.PVALUE" = round(shapiro_pvalue$p.value,3),
        "BARTLETT.PVALUE" = if(exists("bartlett_pvalue")) round(bartlett_pvalue$p.value,3)  else NA,
        "COUNT.CASE"=length(independentVariableData),
        "AVG.CASE" = round(mean(independentVariableData),3),
        "SD.CASE"=round(stats::sd(independentVariableData),3),
        "COUNT.CONTROL"=length(stats::na.omit(dependentVariableData)),
        "AVG.CONTROL"=round(mean(dependentVariableData),3),
        "SD.CONTROL"= round(stats::sd(dependentVariableData),3),
        "RHO"= if(exists("result.cor"))  round(result.cor$estimate) else NA
      )
    } else
    {
      # browser()
      dependentVariableData <- stats::na.omit(tempDataFrame[!is.na(tempDataFrame[,independentVariable]),burdenValue])
      independentVariableData <- stats::na.omit(tempDataFrame[  ,independentVariable])
      local_result <- data.frame (
        "INDIPENDENT.VARIABLE"= independentVariable,
        "ANOMALY" = key$ANOMALY,
        "FIGURE" = key$FIGURE,
        "GROUP" = key$GROUP,
        "SUBGROUP" = key$SUBGROUP,
        "AREA_OF_TEST" = burdenValue,
        "PVALUE" = round(pvalue,3),
        "PVALUEADJ" = pvalueadjusted,
        "TEST" = "SINGLE_AREA",
        "BETA" = if(exists("result.glm")) round(summary( result.glm )$coeff[-1, 1][1],3) else NA,
        "AIC" = if(exists("result.glm")) round(result.glm$aic,3)  else NA,
        "RESIDUALS.SUM" = if(exists("result.glm")) round(sum(result.glm$residuals),3)  else NA,
        "FAMILY" = family_test,
        "transformation" = transformation,
        "COVARIATES" = paste0(covariates,collapse=" "),
        "SHAPIRO.PVALUE" = round(shapiro_pvalue$p.value,3),
        "BARTLETT.PVALUE" = NA,
        "COUNT.CASE"=length(dependentVariableData),
        "AVG.CASE" = round(mean(dependentVariableData),3),
        "SD.CASE"=round(stats::sd(dependentVariableData ),3),
        "COUNT.CONTROL"=length(independentVariableData),
        "AVG.CONTROL"=round(mean(independentVariableData),3),
        "SD.CONTROL"= round(stats::sd(independentVariableData),3),
        "RHO"= if(exists("result.cor"))  round(result.cor$estimate) else NA
      )
    }
    if (exists("result_temp"))
    {
      result_temp <- rbind(result_temp, local_result)
    }
    else
    {
      result_temp <- local_result
    }
  }

  result_temp <- unique(result_temp)


  result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUEADJ"]  <- round(stats::p.adjust(result_temp[result_temp$AREA_OF_TEST=="TOTAL","PVALUE"]  ,method = "BH"),3)
  result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUEADJ"]  <- round(stats::p.adjust(result_temp[result_temp$AREA_OF_TEST!="TOTAL","PVALUE"]  ,method = "BH"),3)

  gc()
  return(result_temp)
}
