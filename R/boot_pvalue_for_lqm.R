bca_pvalue_for_lqm <-
  function(estimate,
           boot_vector,
           model,
           working_data)
  {
    # browser()
    # boot.bca <- coxed::bca(na.omit(boot_vector))
    # p_value <- boot.p_value::boot.p_value(na.omit(boot_vector), type = "bca", theta_null = 0)

    standard_error <- stats::sd(boot_vector, na.rm = T)
    mean_boot_vector <- mean(boot_vector, na.rm = T)
    bias <- estimate - mean_boot_vector
    bias_threshold <- 0.25 * standard_error
    p_value <- 2 * stats::pt(-abs((estimate) / standard_error), (nrow(working_data) - model$edf))
    ci_lower_limit <- estimate + stats::qt(0.05 / 2, (nrow(working_data) - model$edf)) * standard_error
    ci_upper_limit <- estimate + stats::qt(1 - 0.05 / 2, (nrow(working_data) - model$edf)) * standard_error

    #Bias-corrected p-value
    if (abs(bias) > bias_threshold) {
      p_value <- 2 * stats::pt(-abs((2 * estimate - mean_boot_vector) / standard_error), (nrow(working_data) - model$edf))
      ci_lower_limit <-  2 * estimate - mean_boot_vector + stats::qt(0.05 / 2, (nrow(working_data) - model$edf)) * standard_error
      ci_upper_limit <- 2 * estimate - mean_boot_vector + stats::qt(1 - 0.05 / 2, (nrow(working_data) - model$edf)) * standard_error
    }
    # return(c(ci_lower_limit,ci_upper_limit,p_value))

    # DiCiccio, T. J. and B. Efron. (1996). Bootstrap Confidence Intervals. Statistical Science. 11(3):
    # 189–212. https://doi.org/10.1214/ss/1032280214

    # boot.bca <- coxed::bca(na.omit(boot_vector))
    # if(boot.bca[1]<0 & boot.bca[2]>0)
    #   p_value <- 0.06
    # else
    #   p_value <- 0.04
    # return(c(boot.bca,p_value))

  }
