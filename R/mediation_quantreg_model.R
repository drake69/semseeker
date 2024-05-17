# (family_test, tempDataFrame, sig.formula, transformation, plot)
mediation_quantreg_model <- function(family_test,tempDataFrame, sig.formula, transformation, plot) {

  # independent_variable, burden, dependent_variable
  # independent_variable = "age"
  # burden = "methylation"
  # dependent_variable = "t-score","cancer"
  # covariates
  # browser()
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  res <- data.frame()
  tau <- as.numeric(quantreg_params[2])
  res <- data.frame("tau" = tau)
  permutations_test <- as.numeric(quantreg_params[3])
  res$permutations_test <- permutations_test
  permutations <- as.numeric(quantreg_params[4])
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

  mediator_formula <-(paste(parse(text=burden), "~", parse(text=independent_variable)))
  outcome_formula <- (paste(parse(text=dependent_variable), "~", parse(text=independent_variable), "+", parse(text=burden)))
  if(length(covariates)!=0)
  {
    mediator_formula <- (paste(parse(text=burden), "~", parse(text=independent_variable), "+", paste(covariates, collapse = "+")))
    outcome_formula <- (paste(parse(text=dependent_variable), "~", parse(text=independent_variable), "+", parse(text=burden), "+", paste(covariates, collapse = "+")))
  }
  res$mediator_formula <- mediator_formula
  res$outcome_formula <- outcome_formula
  mediator_formula <- as.formula(mediator_formula)
  outcome_formula <- as.formula(outcome_formula)

  # Set the environment of the formulas to the parent frame (global environment)
  # environment(mediator_formula) <- .GlobalEnv
  # environment(outcome_formula) <- .GlobalEnv



  # Fit the quantile regression models at the median (tau quantile)
  m_model <- quantreg::rq(formula = mediator_formula, data = tempDataFrame, tau = tau)
  y_model <- quantreg::rq(formula =  outcome_formula, data = tempDataFrame, tau = tau)

  assign("outcome_formula", outcome_formula, envir = .GlobalEnv)
  assign("mediator_formula", mediator_formula, envir = .GlobalEnv)
  assign("m_model", m_model, envir = .GlobalEnv)
  assign("y_model", y_model, envir = .GlobalEnv)
  assign("burden", burden, envir = .GlobalEnv)
  assign("independent_variable", independent_variable, envir = .GlobalEnv)
  assign("dependent_variable", dependent_variable, envir = .GlobalEnv)
  assign("covariates", covariates, envir = .GlobalEnv)
  assign("permutations_test", permutations_test, envir = .GlobalEnv)
  assign("permutations", permutations, envir = .GlobalEnv)
  assign("tau", tau, envir = .GlobalEnv)


  # # Set the environment of the formulas to the parent frame
  # environment(outcome_formula) <- .GlobalEnv
  # environment(mediator_formula) <- .GlobalEnv
  # environment(tempDataFrame) <- .GlobalEnv
  # environment(m_model) <- .GlobalEnv
  # environment(y_model) <- .GlobalEnv
  # environment(burden) <- .GlobalEnv
  # environment(independent_variable) <- .GlobalEnv
  # environment(dependent_variable) <- .GlobalEnv
  # environment(covariates) <- .GlobalEnv
  # environment(permutations_test) <- .GlobalEnv
  # environment(permutations) <- .GlobalEnv

  if(length(covariates)==0)
  {
    # Perform mediation analysis
    mediation_result <- mediation::mediate(model.m =  eval(m_model), model.y =  eval(y_model), treat = independent_variable,
      mediator = burden, boot = TRUE, sims = permutations_test)

    mediation_summary <- summary(mediation_result)
    if(mediation_summary$d0.p < as.numeric(ssEnv$alpha))
      mediation_result <- mediation::mediate(model.m =  m_model, model.y =  y_model, treat = independent_variable,
        mediator = burden, boot = TRUE, sims = permutations)
  }
  else
  {
    # eval(bquote(mediate(model.m = mediator_model, model.y = outcome_model,
    #   treat = .(treatment), mediator = .(mediator),
    #   boot = TRUE, sims = 1000)))

    # Perform mediation analysis
    mediation_result <- mediation::mediate(model.m =  m_model, model.y =  y_model, treat = independent_variable,
      mediator = burden, boot = TRUE, sims = permutations_test,
      covariates = tempDataFrame[ ,covariates])

    mediation_summary <- summary(mediation_result)
    if(mediation_summary$d0.p < as.numeric(ssEnv$alpha))
      mediation_result <- mediation::mediate(model.m =  m_model, model.y =  y_model, treat = independent_variable,
        mediator = burden, boot = TRUE,
        sims = permutations,
        covariates = tempDataFrame[ ,covariates])
  }

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

  # remove assigned variables from GlobalEnv
  rm(outcome_formula, mediator_formula, m_model, y_model, burden, independent_variable, dependent_variable, covariates,
    permutations_test, permutations, tempDataFrame, mediation_result, mediation_summary, tau)

  # Return the summary of the mediation analysis
  summary(res)
}

