# (family_test, tempDataFrame, sig.formula, transformation, plot)
mediation_linear_model <- function(family_test,tempDataFrame, sig.formula, transformation, plot, samples_sql_condition=samples_sql_condition, area, subarea) {

  model_params <- unlist(strsplit(as.character(family_test),"_"))
  permutations_test <- as.numeric(model_params[2])
  res <- data.frame(permutations_test = permutations_test)
  permutations <- as.numeric(model_params[3])
  res$permutations <- permutations

  # decompose the formula, to move the mediator
  # the expected formula is "mediator ~ outcome + treatment + covariates"
  vars <- sig.formula_vars(sig.formula)
  mediator <- vars[[1]] # burdnValue
  outcome <- vars[[2]] # dependent variable
  treatment <- vars[[3]][1]

  res$mediator <- mediator
  res$outcome <- outcome
  res$treatment <- treatment

  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- sapply(covariates, as.character)
  } else
    covariates <- c()

  df_mediator <- data.frame(tempDataFrame)
  df_mediator$mediator <- df_mediator[,mediator]
  df_mediator$treatment <- df_mediator[,treatment]
  if(length(covariates) != 0)
    df_mediator <- df_mediator[,c("treatment", "mediator",covariates)]
  else
    df_mediator <- df_mediator[,c("treatment", "mediator")]

  df_outcome <- data.frame(data)
  df_outcome$outcome <- df_outcome[,outcome]
  df_outcome$mediator <- df_outcome[,mediator]
  df_outcome$treatment <- df_outcome[,treatment]
  if(length(covariates) != 0)
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment",covariates)]
  else
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment")]

  perms <- sort(unique(c(permutations_test, permutations)))

  for (i in 1:length(perms)) {

    env <- new.env()
    env$mediate_env <- mediation::mediate
    perm <- perms[i]
    mediation_result <- eval(expr = quote(mediate_env(
      model.m = lm(as.formula(paste("mediator", "~", ".")), data = df_mediator),
      model.y = lm(as.formula(paste("outcome", "~", ".")), data = df_outcome),
      treat = "treatment",
      mediator = "mediator",
      boot = TRUE,
      sims = perm
    )), envir = env)
    res$pvalue <- max(mediation_result[1,4], mediation_result[2,4], mediation_result[3,4], mediation_result[4,4])
    res$permutations <- perm
    # Extract key components from the summary
    res$acme <- mediation_result[1,1]
    res$acme_ci_start <- mediation_result[1,2]
    res$acme_ci_end <- mediation_result[1,3]
    res$acme_pvalue <- mediation_result[1,4]

    res$ade <- mediation_result[2,1]
    res$ade_ci_start <- mediation_result[2,2]
    res$ade_ci_end <- mediation_result[2,3]
    res$ade_pvalue <- mediation_result[2,4]

    res$total_effect <- mediation_result[3,1]
    res$total_effect_ci_start <- mediation_result[3,2]
    res$total_effect_ci_end <- mediation_result[3,3]
    res$total_effect_pvalue <- mediation_result[3,4]

    res$prop_mediated <- mediation_result[4,1]
    res$prop_mediated_ci_start <- mediation_result[4,2]
    res$prop_mediated_ci_end <- mediation_result[4,3]
    res$r_model <- family_test
    if(res$pvalue < as.numeric(ssEnv$alpha))
      next
  }

  # remove environment
  rm(env)


  # Return the summary of the mediation analysis
  return(res)
}

