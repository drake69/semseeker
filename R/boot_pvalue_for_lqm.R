bca_pvalue_for_lqm <-
  function(estimate,
           boot_vector,
           model,
           n_elements)
  {
    # browser()
    # boot.bca <- coxed::bca(na.omit(boot_vector))
    # p_value <- boot.p_value::boot.p_value(na.omit(boot_vector), type = "bca", theta_null = 0)

    if(stats::shapiro.test(na.omit(boot_vector))$p.value > 0.05)
    {
      #boot vector is a normal distribution
      standard_error <- stats::sd(boot_vector, na.rm = T)
      mean_boot_vector <- mean(boot_vector, na.rm = T)
      bias <- estimate - mean_boot_vector
      bias_threshold <- 0.25 * standard_error
      p_value <- 2 * stats::pt(-abs((estimate) / standard_error), (n_elements - model$edf))
      ci_lower_limit <- estimate + stats::qt(0.05 / 2, (n_elements - model$edf)) * standard_error
      ci_upper_limit <- estimate + stats::qt(1 - 0.05 / 2, (n_elements - model$edf)) * standard_error

      #Bias-corrected p-value
      if (abs(bias) > bias_threshold) {
        p_value <- 2 * stats::pt(-abs((2 * estimate - mean_boot_vector) / standard_error), (n_elements - model$edf))
        ci_lower_limit <-  2 * estimate - mean_boot_vector + stats::qt(0.05 / 2, (n_elements - model$edf)) * standard_error
        ci_upper_limit <- 2 * estimate - mean_boot_vector + stats::qt(1 - 0.05 / 2, (n_elements - model$edf)) * standard_error
      }
      return(c(ci_lower_limit,ci_upper_limit,p_value))
    }
    else
    {
      # DiCiccio, T. J. and boot_vector. Efron. (1996). Bootstrap Confidence Intervals. Statistical Science. 11(3):
      # 189–212. https://doi.org/10.1214/ss/1032280214
      boot.bca <- coxed::bca(na.omit(boot_vector))
      if(boot.bca[1] <0 & boot.bca[2]>0)
        p_value <-1
      else
        p_value<-0
      return(c(boot.bca,p_value))
    }

  }




# BCApval<-function(boot_vector, estimate, working_data){
#
#   alpha_seq <- seq(1e-16, 1-1e-16, 5e-6)
#   #BCA method
#   #Desired quantiles
#
#   u <- cbind(alpha_seq/2, 1-alpha_seq/2)
#
#   #Compute constants
#   z0 <- qnorm(mean(boot_vector <= estimate))
#   zu <- qnorm(u)
#
#   group_all <- working_data$soggetto
#   group_unique <- unique(group_all)
#
#   n<-length(group_unique)
#   I <- rep(NA, n)
#   i<-1
#   for(i in 1:n){
#     #Remove ith working_data point
#     xnew <- working_data[-which(working_data$soggetto==group_unique[i]),]
#
#     #Estimate beta dal modello
#     rqrja <- quantreg(modello)
#
#     I[i] <- (n-1)*(estimate -  rqrja$beta)
#   }
#   #Estimate a
#   a <- (sum(I^3)/sum(I^2)^1.5)/6
#
#   #Adjusted quantiles
#   u_adjusted <- pnorm(z0 + (z0+zu)/(1-a*(z0+zu)))
#
#   #Accelerated Bootstrap CI
#   Bca1<-quantile(boot_vector, u_adjusted[,1])
#   Bca2<-quantile(boot_vector, u_adjusted[,2])
#   Bca<-cbind(Bca1,Bca2)
#   pval <- alpha_seq[which.min(0 >= Bca[,1] & 0 <= Bca[,2])]
#   return(pval)
# }
