apply_stat_model_sig_formula <- function (family_test, burdenValue, independent_variable, covariates)
{
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

  return (sig.formula)
}
