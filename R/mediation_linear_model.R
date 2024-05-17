# (family_test, tempDataFrame, sig.formula, transformation, plot)
mediation_linear_model <- function(family_test,tempDataFrame, sig.formula, transformation, plot) {

  # independent_variable, burden, dependent_variable
  # independent_variable = "age"
  # burden = "methylation"
  # dependent_variable = "t-score","cancer"
  # covariates

  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  permutations_test <- as.numeric(quantreg_params[2])
  res <- data.frame(permutations_test = permutations_test)
  permutations <- as.numeric(quantreg_params[3])
  res$permutations <- permutations

  # decompose the formula, to move the mediator
  # the expected formula is "burden ~ dependent_variable + independent_variable + covariates"
  vars <- sig.formula_vars(sig.formula)
  burden <- vars[[1]]
  dependent_variable <- vars[[2]]
  independent_variable <- vars[[3]][1]
  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- sapply(covariates, as.character)
  } else
  {
    covariates <- c()
  }

  if(length(covariates)==0)
  {
    mediator_formula <- as.formula(paste(parse(text=burden), "~", parse(text=independent_variable)))
    outcome_formula <- as.formula(paste(parse(text=dependent_variable), "~", parse(text=independent_variable), "+", parse(text=burden)))
  }
  else
  {
    mediator_formula <- as.formula(paste(parse(text=burden), "~", parse(text=independent_variable), "+", paste(covariates, collapse = "+")))
    outcome_formula <- as.formula(paste(parse(text=dependent_variable), "~", parse(text=independent_variable), "+", parse(text=burden), "+", paste(covariates, collapse = "+")))
  }

  # Fit the quantile regression models at the median (tau quantile)
  m_model <- lm(formula = mediator_formula, data = tempDataFrame)
  y_model <- lm(formula =  outcome_formula, data = tempDataFrame)

  # Perform mediation analysis
  mediation_result <- mediation::mediate(model.m =  m_model, model.y =  y_model, treat = independent_variable, mediator = burden, boot = TRUE, sims = permutations_test,
    covariates = tempDataFrame[ ,covariates])

  if(mediation_summary$d0.p < ssEnv$alpha)
    mediation_result <- mediation::mediate(model.m =  m_model, model.y =  y_model, treat = independent_variable, mediator = burden, boot = TRUE, sims = permutations,
      covariates = tempDataFrame[ ,covariates])

  # Return the summary of the mediation analysis
  mediation_summary <- summary(mediation_result)

  # Extract key components from the summary
  res$acme <- mediation_summary$d0
  res$acme_ci_start <- mediation_summary$d0.ci[1]
  res$acme_ci_end <- mediation_summary$d0.ci[2]
  res$acme_pvalue <- mediation_summary$d0.p

  res$ade <- mediation_summary$z0
  res$ade_ci_start <- mediation_summary$z0.ci[1]
  res$ade_ci_end <- mediation_summary$z0.ci[2]
  res$ade_pvalue <- mediation_summary$z0.p

  res$total_effect <- mediation_summary$tau.coef
  res$total_effect_ci_start <- mediation_summary$tau.ci[1]
  res$total_effect_ci_end <- mediation_summary$tau.ci[2]
  res$total_effect_pvalue <- mediation_summary$tau.p

  res$prop_mediated <- mediation_summary$n0
  res$prop_mediated_ci_start <- mediation_summary$n0.ci[1]
  res$prop_mediated_ci_end <- mediation_summary$n0.ci[2]

  # Return the summary of the mediation analysis
  summary(res)
}

