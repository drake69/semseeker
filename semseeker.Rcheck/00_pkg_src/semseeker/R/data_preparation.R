#' Title
#'
#' @param family_test test or regression to apply
#' @param transformation transformation to apply to data
#' @param tempDataFrame data frame to use for test/regression
#' @param independent_variable regressor
#' @param g_start starting column of the dataframe
#' @param dototal boolean to calculate the total burden test/regression
#' @param covariates vector of covariates to be found in the sample sheet
#' @param depth_analysis 1 only sample, 2 chr, 3 alle genomic areas
#'
data_preparation <- function(family_test,transformation,tempDataFrame, independent_variable, g_start, dototal, covariates, depth_analysis)
{

  ssEnv <- get_session_info()

  transformation <- as.character(transformation)
  originalDataFrame <- tempDataFrame
  tempDataFrame <- as.data.frame(sapply(tempDataFrame, as.numeric))

  independent_variable1stLevel <- NA
  independent_variable2ndLevel <- NA
  test_factor <- as.factor(tempDataFrame[, independent_variable])
  if(length(levels(test_factor))<4)
  {
    tempDataFrame[, independent_variable] <- as.factor(tempDataFrame[, independent_variable])
    }
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

  if(family_test == "log")
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
    suppressWarnings(
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
    stop("ERROR: I'm stopping here data are not the same size, file a bug!")

  colnames(tempDataFrame) <- df_colnames
  # after the transformation some data could be missed
  lost_cols <- colSums(apply(tempDataFrame,2,is.nan))!=0
  lostDataFrame <-  colnames(tempDataFrame)[lost_cols]
  if(sum(lost_cols)!=0)
    utils::write.csv2(lostDataFrame, file.path(ssEnv$session_folder,paste("lost_data_",transformation,"_",stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),".log", sep="")))

  #  we want to preserve the NA in the independent variables to be removed by the models
  tempDataFrame[apply(tempDataFrame,2,is.nan)] <- 0

  result <- list(tempDataFrame, independent_variable1stLevel, independent_variable2ndLevel)
  names(result) <- c("tempDataFrame", "independent_variable1stLevel","independent_variable2ndLevel")
  return (result)
}
