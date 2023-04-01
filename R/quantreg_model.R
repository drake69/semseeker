#' Title
#'
#' @param sig.formula
#' @param df
#' @param tau
#' @param lqm_control
#'
#' @return
#' @export
#'
#' @examples
compute_quantreg_beta_boot_np <- function(sig.formula,df, tau, lqm_control)
{
  cols <- colnames(df)
  tempDataFrame <- Rfast::colShuffle(as.matrix(df))
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
  beta_value <- summary(model)$tTable[2,"Value"]
  return(beta_value)
}

#' Title
#'
#' @param sig.formula
#' @param tau
#' @param localDataFrame
#'
#' @return
#' @export
#'
#' @examples
compute_qr_beta_boot_p <- function(sig.formula, tau, localDataFrame) {
  fit <- quantreg::rq(formula =  sig.formula,data = as.data.frame(localDataFrame),  tau = tau)
  coef <-as.data.frame(summary(fit, se = "boot")$coefficients)[2,"Value"]
  pval <- as.data.frame(summary(fit, se = "boot")$coefficients)[2, "Pr(>|t|)"]
  return(list(beta = coef, pval = pval))
}


#' Title
#'
#' @param family_test
#' @param sig.formula
#' @param tempDataFrame
#' @param independent_variable
#' @param boot_success
#' @param tests_count
#'
#' @return
#' @export
#'
#' @examples
quantreg_model <- function(family_test, sig.formula, tempDataFrame, independent_variable, boot_success, tests_count)
{
  n_permutations <- NA
  ci.lower.adjusted <- NA
  ci.upper.adjusted <- NA

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  # quantreg_params template quantreg + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  if(length(quantreg_params)<5)
  {
    # just apply quantile regression only 2 params quantreg + quantile
    if(length(quantreg_params)<4)
    {
      tau = as.numeric(quantreg_params[2])
      model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
      pvalue <- summary(model)$tTable[2,"Pr(>|t|)"]
      std.error <- summary(model)$tTable[2,"Std. Error"]
      beta_value <- summary(model)$tTable[2,"Value"]
      ci.lower <- summary(model)$tTable[2,"lower bound"]
      ci.upper <- summary(model)$tTable[2,"upper bound"]
      r_model <- "lqmm_lqm"
    }
    else
    {
      # apply bootstrap to obtain beta using quantreg
      # quantreg + quantile + first_round_of_permutations + second_round_of_permutations
      # Define function to compute p-value and beta regression coefficient
      tau = as.numeric(quantreg_params[2])
      n_permutations_test <- as.numeric(quantreg_params[3])
      # Compute beta and p-value for n_permutations replications
      results <- replicate(n_permutations_test, compute_qr_beta_boot_p(sig.formula, tau, localDataFrame=tempDataFrame))
      # Compute average beta and p-value
      pvalue <- max(unlist(t(results)[,"pval"]))
      n_permutations <- as.numeric(quantreg_params[4])
      if(pvalue < 0.05 && n_permutations_test < n_permutations)
      {
        results <- replicate(n_permutations, compute_qr_beta_boot_p(sig.formula, tau, localDataFrame=tempDataFrame))
        # Compute average beta and p-value
        pvalue <- max(unlist(t(results)[,"pval"]))
      }
      beta_value <- mean(unlist(t(results)[,"beta"]))
      r_model <- "quantreg_rq"
    }
  }
  if(length(quantreg_params)==5)
  {
    # use lqmm package and apply quantreg with bootstrap and confidence interval of regression beta
    # quantreg + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
    n_permutations_test <- as.numeric(quantreg_params[3])
    n_permutations <- as.numeric(quantreg_params[4])
    tau <- as.numeric(quantreg_params[2])
    conf.level <- as.numeric(quantreg_params[5])
    model.x <-  suppressMessages(lqmm::lqm( formula = sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control))
    if(n_permutations > n_permutations_test)
    {
      model.x.boot <- suppressMessages(lqmm::boot(model.x, R = n_permutations_test))
      beta_value <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
      std.error <- summary(model.x.boot)$tTable[2,"Std. Error"]
      tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
      colnames(tt) <- colnames(model.x.boot)
      boot_vector <- stats::na.omit(tt[,independent_variable])
      boot.bca <- quantreg_summary(boot_vector, beta_value, conf.level = conf.level)
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
      std.error <- summary(model.x.boot)$tTable[2,"Std. Error"]
      tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
      colnames(tt) <- colnames(model.x.boot)
      boot_vector <- stats::na.omit(tt[,independent_variable])
      boot.bca <- quantreg_summary(boot_vector, beta_value, conf.level = conf.level)
      boot.bca.adjusted <- quantreg_summary(boot_vector, beta_value,  boot_success = boot_success, tests_count=tests_count, conf.level = conf.level)
      ci.lower.adjusted <-  boot.bca.adjusted[1]
      ci.upper.adjusted <- boot.bca.adjusted[2]
    }

    ci.lower <- boot.bca[1]
    ci.upper <- boot.bca[2]
    pvalue <- boot.bca[3]
    r_model <- "lqmm_lqm"
  }

  if(length(quantreg_params)==6)
  {
    tau <- as.numeric(quantreg_params[2])
    n_permutations_test <- as.numeric(quantreg_params[3])
    n_permutations <- as.numeric(quantreg_params[4])
    conf.level <- as.numeric(quantreg_params[5])
    parametric <- as.numeric(quantreg_params[6]=="p")

    model <- lqmm::lqm(sig.formula, tau =tau, data = as.data.frame(tempDataFrame), na.action = stats::na.omit, control = lqm_control)
    beta_value <- summary(model)$tTable[2,"Value"]
    std.error <- summary(model)$tTable[2,"Std. Error"]
    ci.lower <- summary(model)$tTable[2,"lower bound"]
    ci.upper <- summary(model)$tTable[2,"upper bound"]
    pvalue <- summary(model)$tTable[2,"Pr(>|t|)"]

    if(n_permutations > n_permutations_test)
    {
      boot_beta <- replicate(n_permutations_test, compute_quantreg_beta_boot_np(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control))
      boot.bca <- quantreg_summary(boot_beta, beta_value, conf.level = conf.level)
      ci.lower <- boot.bca[1]
      ci.upper <- boot.bca[2]
      pvalue <- boot.bca[3]
    }
    if(pvalue > 0.05)
    {
      n_permutations <- n_permutations_test
    }
    else
    {
      boot_beta <- replicate(n_permutations, compute_quantreg_beta_boot_np(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control))
      boot.bca <- quantreg_summary(boot_beta, beta_value, conf.level = conf.level)
      ci.lower <- boot.bca[1]
      ci.upper <- boot.bca[2]
      pvalue <- boot.bca[3]
    }
    r_model <- "lqmm_lqm"
  }

  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA


  return (data.frame(ci.lower,ci.upper, pvalue, beta_value,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations,ci.lower.adjusted,ci.upper.adjusted))

}
