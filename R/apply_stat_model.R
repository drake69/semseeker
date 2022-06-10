#' @importFrom doRNG %dorng%
apply_stat_model <- function(tempDataFrame, g_start, family_test, covariates = NULL, key, transformation, dototal, logFolder, independent_variable, depth_analysis=3)
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


  independent_variables <- tempDataFrame[,1:(g_start-1)]
  burden_values <- as.data.frame(sapply(tempDataFrame[,g_start:ncol(tempDataFrame)], as.numeric))

  if( ncol(tempDataFrame) - g_start  > 1) {
    burden_values <- burden_values[, colSums(burden_values) > 0]
  }

  df_colnames <- colnames(tempDataFrame)
  if( !is.null(dim(burden_values))  & dototal) {
    # browser()
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

  # #browser()
  if(family_test != "poisson")
    burden_values <- burden_values + 0.001

  transformation <- as.character(transformation)
  if(is.null(transformation) | length(transformation)==0 | is.na(transformation))
    transformation <- "none"

  df_values_orig <- burden_values
  try(
    {
      burden_values = switch(
        as.character(transformation),
        "scale" = as.data.frame(apply(burden_values,2,scale)),
        "log" = log(burden_values),
        "log2" = log2(burden_values),
        "log10"= log10(burden_values),
        "exp" = exp(burden_values),
        "johnson" =  Johnson::RE.Johnson(burden_values)$transformed,
        "none" = burden_values,
        burden_values
      )
    }
  )

  if(grepl("quantile", transformation))
  {
    qq <- unsplit(strsplit(as.character(transformation),"_"))[2]
    df_values_temp <- as.data.frame(apply( burden_values,2,function(x) ggplot2::cut_number(x, n=qq)))
    colnames(df_values_temp) <- colnames(burden_values)
  }

  if(setequal(burden_values,df_values_orig) & transformation !="none")
    transformation <- paste0("NA_", transformation, sep="")


  if(family_test!="binomial" & family_test!="wilcoxon" & family_test!="t.test" & family_test!="poisson")
  {
    independent_variableValues <- as.data.frame(as.numeric(tempDataFrame[, independent_variable]))
    independent_variableValuesOrig <- as.data.frame(as.numeric(tempDataFrame[, independent_variable]))
    try(
      {
        independent_variableValues = switch(
          as.character(transformation),
          "scale" = as.data.frame(apply(independent_variableValues,2,scale)),
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

    if(setequal(burden_values,df_values_orig) & transformation !="none")
      transformation <- paste0("NA_", transformation, sep="")
    else
    {
      tempDataFrame[, independent_variable] <- independent_variableValues
      }
  }

  tempDataFrame <- data.frame(independent_variables, burden_values)
  # browser()
  colnames(tempDataFrame) <- df_colnames
  cols <- colnames(tempDataFrame)
  iters <- length(cols)

  g <- 0
  to_export <- c("cols", "family_test", "covariates", "independent_variable", "tempDataFrame",
                 "independent_variable1stLevel", "independent_variable2ndLevel",
                 "key", "transformation","BCApval")
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
      bartlett_pvalue <- stats::bartlett.test( stats::as.formula(paste0(burdenValue,"~", independent_variable, sep="")),
                                               data= as.data.frame(tempDataFrame) )
    }

    if(family_test=="gaussian" | family_test=="spearman" | family_test=="kendall" | family_test=="pearson")
    {
      # browser()
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
      # browser()
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- independent_variable
      } else
      {
        covariates_model <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
      }
      sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
      lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
      quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
      n_permutations_test <- as.numeric(quantreg_params[3])
      n_permutations <- as.numeric(quantreg_params[4])
      tau <- as.numeric(quantreg_params[2])

      model.x <-  lqmm::lqm(sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control)
      if(n_permutations > n_permutations_test)
      {
        model.x.boot <- lqmm::boot(model.x, R = n_permutations_test)
        beta_full <- summary(model.x.boot)[independent_variable,"Value"]
        tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
        colnames(tt) <- colnames(model.x.boot)
        boot_vector <- stats::na.omit(tt[,independent_variable])
        boot.bca <- BCApval(boot_vector, beta_full, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable)
      }

      if(boot.bca[1]<0 & boot.bca[2]>0)
      {
        n_permutations <- n_permutations_test
      }
      else
      {
        model.x.boot <- lqmm::boot(model.x, R = n_permutations)
        beta_full <- summary(model.x.boot)[independent_variable,"Value"]
        tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
        colnames(tt) <- colnames(model.x.boot)
        boot_vector <- stats::na.omit(tt[,independent_variable])
        boot.bca <- BCApval(boot_vector, beta_full, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable)
      }

      ci.lower <- boot.bca[1]
      ci.upper <- boot.bca[2]
      pvalue <- boot.bca[3]
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
    # browser()
    # beta <- exp(summary( result_glm )$coeff[-1, 1][1])
    # result_temp <-

    if(family_test!="gaussian" & family_test!="spearman" & family_test!="pearson" & family_test!="kendall" & !grepl("quantreg",family_test))
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
        "RHO"= if(exists("result_cor"))  round(result_cor$estimate,4) else NA,
        "CI.LOWER"= if(exists("ci.lower")) ci.lower else NA,
        "CI.UPPER"= if(exists("ci.upper")) ci.upper else NA,
        "N.PERMUTATIONS" =  NA
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
        "RHO"= if(exists("result_cor"))  round(result_cor$estimate,4) else NA,
        "CI.LOWER"= if(exists("ci.lower")) ci.lower else NA,
        "CI.UPPER"= if(exists("ci.upper")) ci.upper else NA,
        "N.PERMUTATIONS" = if(exists("n_permutations")) n_permutations else NA
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
