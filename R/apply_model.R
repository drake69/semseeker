
apply_model <- function(df, g_start, family, covariates, key, transformation, dototal, logFolder)
{

  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
  outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
  # print(outFile)
  computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computation_cluster)

  # options(digits = 22)
  parallel::clusterExport(envir=environment(), cl = computation_cluster, varlist =c())

  df$Sample_Group[df$Sample_Group=="Control"] <- 0
  df$Sample_Group[df$Sample_Group=="Case"] <- 1
  df$Sample_Group <- as.numeric(df$Sample_Group)

  df_head <- df[,1:(g_start-1)]
  df_values <- sapply(df[,g_start:ncol(df)], as.numeric)

  if( ncol(df) - g_start  > 1) {
    df_values <- df_values[, colSums(df_values) > 0]
  }

  if( !is.null(dim(df_values))  & dototal) {
    sum_area <- apply(df_values, 1, sum)
    df_values <- cbind(df_values, data.frame("TOTAL"=sum_area))
  }

  if(family != "poisson")
    df_values <- df_values + 0.001

  transformation_label <- "NONE"
  if(!is.null(transformation))
  {
    transformation_label <- as.character(substitute(transformation))
    df_values <- transformation(df_values)
  }

  df <- cbind(df_head, df_values)
  cols <- colnames(df)
  iters <- length(cols)

  result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind) %dopar%
  # for (g in g_start:iters)
  {
    #g <- 2

    if (family=="poisson")
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- cols[g]
      } else
      {
        covariates_model <- paste(paste0(c("Sample_Group", covariates),collapse="+", sep=""))
      }
      sig.formula <- as.formula(paste0(cols[g],"~", covariates_model, sep=""))
    } else
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        covariates_model <- cols[g]
      } else
      {
        covariates_model <- paste(paste0(c(cols[g], covariates),collapse="+", sep=""))
      }
      sig.formula <- as.formula(paste0("Sample_Group","~", covariates_model, sep=""))
    }

    result.glm  <- glm( sig.formula, family = family, data = df)
    pvalue <- summary( result.glm )$coeff[-1, 4][1]
    pvalueadjusted <- pvalue
    # beta <- exp(summary( result.glm )$coeff[-1, 1][1])
    # result_temp <-
      data.frame (
        "ANOMALY" = key$ANOMALY,
        "FIGURE" = key$FIGURE,
        "GROUP" = key$GROUP,
        "SUBGROUP" = key$SUBGROUP,
        "AREA_OF_TEST" = cols[g],
        "PVALUE" = pvalue,
        "PVALUEADJ" = pvalueadjusted,
        "TEST" = "SINGLE_AREA",
        "BETA" = summary( result.glm )$coeff[-1, 1][1],
        "AIC" = result.glm$aic,
        "RESIDUALS.SUM" = sum(result.glm$residuals),
        "FAMILY" = family,
        "TRANSFORMATION" = transformation_label
      )
  }
  result_temp$PVALUEADJ <- p.adjust(result_temp$PVALUE,method = "BH")
  result_temp$TRANSFORMATION <- transformation_label
  return(result_temp)
}
