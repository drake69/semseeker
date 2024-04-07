#' Title
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param tau tau at which apply the wuantile regression
#' @param lqm_control specification of the lqmm package
#'
compute_quantreg_signal_boot_np <- function(sig.formula,df, tau, lqm_control)
{
  # browser()
  cols <- colnames(df)
  # tempDataFrame <- Rfast::colShuffle(as.matrix(df))
  indepVar <- as.character(all.vars(sig.formula)[-1][1])
  tempDataFrame <- df
  tempDataFrame[ ,indepVar] <- sample(tempDataFrame[,indepVar])
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  suppressMessages({
    model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
  })
  summary_qr <- suppressMessages(summary(model)$tTable)
  signal_value <- summary_qr[2,"Value"]
  return(signal_value)
}

#' Title
#'
#' @param sig.formula formula to use for regression application
#' @param tau tau to apply the quantile regression
#' @param localDataFrame dataframe to apply th regression model
#'
compute_qr_signal_boot_p <- function(sig.formula, tau, localDataFrame) {
  suppressMessages({
    fit <- quantreg::rq(formula =  sig.formula,data = as.data.frame(localDataFrame),  tau = tau)
  })
  coef <-as.data.frame(summary(fit, se = "boot")$coefficients)[2,"Value"]
  pval <- as.data.frame(summary(fit, se = "boot")$coefficients)[2, "Pr(>|t|)"]
  return(list(signal = coef, pval = pval))
}


#' Title
#'
#' @param family_test family lqmm, quantreg
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param boot_success number of success tests to calculate corrected confidence interval
#' @param tests_count count of total executed tests
#'
quantreg_model <- function(family_test, sig.formula, tempDataFrame, independent_variable, boot_success, tests_count)
{
  n_permutations <- NA
  ci.lower.adjusted <- NA
  ci.upper.adjusted <- NA
  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  signal_value <- NA
  pvalue <- NA
  boot.bca <- NA

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 10000, verbose = F )
  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  # quantreg_params template quantreg + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  if(length(quantreg_params)<5)
  {
    # just apply quantile regression only 2 params quantreg + quantile
    if(length(quantreg_params)<4)
    {
      tau = as.numeric(quantreg_params[2])
      try({
        model <- lqmm::lqm(sig.formula, tau =tau, data = tempDataFrame, na.action = stats::na.omit, control = lqm_control)
        summary_qr <- suppressMessages(summary(model)$tTable)
        pvalue <- summary_qr[2,"Pr(>|t|)"]
        std.error <- summary_qr[2,"Std. Error"]
        signal_value <- summary_qr[2,"Value"]
        ci.lower <- summary_qr[2,"lower bound"]
        ci.upper <- summary_qr[2,"upper bound"]
      })
      r_model <- "lqmm_lqm"
    }
    else
    {
      # apply bootstrap to obtain signal using quantreg
      # quantreg + quantile + first_round_of_permutations + second_round_of_permutations
      # Define function to compute p-value and signal regression coefficient
      tau = as.numeric(quantreg_params[2])
      n_permutations_test <- as.numeric(quantreg_params[3])
      # Compute signal and p-value for n_permutations replications
      results <- replicate(n_permutations_test, compute_qr_signal_boot_p(sig.formula, tau, localDataFrame=tempDataFrame))
      # Compute average signal and p-value
      pvalue <- max(unlist(t(results)[,"pval"]))
      n_permutations <- as.numeric(quantreg_params[4])
      if(pvalue < 0.05 && n_permutations_test < n_permutations)
      {
        results <- replicate(n_permutations, compute_qr_signal_boot_p(sig.formula, tau, localDataFrame=tempDataFrame))
        # Compute average signal and p-value
        pvalue <- max(unlist(t(results)[,"pval"]))
      }
      signal_value <- mean(unlist(t(results)[,"signal"]))
      r_model <- "quantreg_rq"
    }
  }
  if(length(quantreg_params)==5)
  {
    # use lqmm package and apply quantreg with bootstrap and confidence interval of regression signal
    # quantreg + quantile + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
    n_permutations_test <- as.numeric(quantreg_params[3])
    n_permutations <- as.numeric(quantreg_params[4])
    tau <- as.numeric(quantreg_params[2])
    conf.level <- as.numeric(quantreg_params[5])
    suppressMessages({
      model.x <-  suppressMessages(lqmm::lqm( formula = sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control))
    })
    if(n_permutations > n_permutations_test)
    {
      model.x.boot <- suppressMessages(lqmm::boot(model.x, R = n_permutations_test))
      signal_value <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
      std.error <- summary(model.x.boot)[2,"Std. Error"]
      tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
      colnames(tt) <- colnames(model.x.boot)
      boot_vector <- stats::na.omit(tt[,independent_variable])
      boot.bca <- quantreg_summary(boot_vector, signal_value, conf.level = conf.level)
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
      signal_value <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
      std.error <- summary(model.x.boot)[2,"Std. Error"]
      tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
      colnames(tt) <- colnames(model.x.boot)
      boot_vector <- stats::na.omit(tt[,independent_variable])
      boot.bca <- quantreg_summary(boot_vector, signal_value, conf.level = conf.level)
      boot.bca.adjusted <- quantreg_summary(boot_vector, signal_value,  boot_success = boot_success, tests_count=tests_count, conf.level = conf.level)
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
    suppressMessages({
      model <- lqmm::lqm(sig.formula, tau =tau, data = as.data.frame(tempDataFrame), na.action = stats::na.omit, control = lqm_control)
    })
    summary_qr <- suppressMessages(summary(model)$tTable)
    signal_value <- summary_qr[2,"Value"]
    std.error <- summary_qr[2,"Std. Error"]
    ci.lower <- summary_qr[2,"lower bound"]
    ci.upper <- summary_qr[2,"upper bound"]
    pvalue <- summary_qr[2,"Pr(>|t|)"]

    if(n_permutations > n_permutations_test)
    {
      boot_beta <- suppressMessages(replicate(n_permutations_test, compute_quantreg_signal_boot_np(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control)))
      boot.bca <- quantreg_summary(boot_beta, signal_value, conf.level = conf.level)
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
      boot_beta <- suppressMessages(replicate(n_permutations, compute_quantreg_signal_boot_np(sig.formula,as.data.frame(tempDataFrame), tau, lqm_control)))
      boot.bca <- quantreg_summary(boot_beta, signal_value, conf.level = conf.level)
      ci.lower <- boot.bca[1]
      ci.upper <- boot.bca[2]
      pvalue <- boot.bca[3]
    }
    r_model <- "lqmm_lqm"
  }



  return (data.frame(ci.lower,ci.upper, pvalue, signal_value,aic_value,residuals,shapiro_pvalue,r_model,std.error,n_permutations,ci.lower.adjusted,ci.upper.adjusted))

}
