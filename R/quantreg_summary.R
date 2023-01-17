#' Quantile regression result value, confidence interval and pvalue
#'
#' @param boot_vector vector of boot statistc beta regression
#' @param estimate beta regression
#' @param working_data data to regress
#' @param sig.formula formula for model
#' @param tau quantile to regress at
#' @param independent_variable name of indenpendent variable
#' @param lqm_control controls of lqmm packages
#' @param conf.level confidence intervals alpha level
#' @param boot_success number of success respecting the null hypothesis
#' @param tests_count how many tests were done
#'
#' @return ci and pvalue with BCA method
#' @importFrom doRNG %dorng%
quantreg_summary <-function(boot_vector, estimate, working_data, sig.formula, tau, independent_variable, lqm_control, conf.level , boot_success = 0, tests_count=1 ){


  # #BCA method
  # #Desired quantiles
  # alpha<-0.05
  # alpha_seq <- seq(1e-16, 1-1e-16, 5e-6)
  #
  # alpha_interval <- c(alpha/2, 1-alpha/2)
  # alpha_interval_seq <- cbind(alpha_seq/2, 1-alpha_seq/2)
  #
  # #Compute constants
  # inverse_cumulative_normal_standard_alpha_of_mean_difference <- stats::qnorm(mean(boot_vector <= estimate), mean = 0, sd = 1)
  # inverse_cumulative_normal_standard_alpha <- stats::qnorm(alpha_interval)
  # inverse_cumulative_normal_standard_alpha_seq <- stats::qnorm(alpha_interval_seq)
  #
  # n<-length(working_data)
  # I <- rep(NA, n)
  # # jack <- I
  # i<- 0
  # to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable")
  # # I <- foreach::foreach(i = 1:n, .combine = c, .export = to_export) %dorng%
  # for(i in 1:n)
  # {
  #   #Remove ith working_data point
  #   working_data_new <- working_data[-i,]
  #
  #   #Estimate beta dal modello
  #   fit.lqm <- suppressMessages(lqmm::lqm( formula = sig.formula, data = working_data_new, tau = tau, control = lqm_control, fit = TRUE))
  #   fit.res <- suppressMessages(summary(fit.lqm)$tTable)
  #
  #   # jack[i] <-  fit.res[independent_variable,"Value"]
  #   I[i]<-(n-1)*(estimate -  fit.res[independent_variable,"Value"] )
  # }
  #
  # #Estimate acceleration
  # acceleration <- (sum(I^3)/sum(I^2)^1.5)/6
  #
  #
  # #Adjusted quantiles
  # u_adjusted_seq <- stats::pnorm(inverse_cumulative_normal_standard_alpha_of_mean_difference + (inverse_cumulative_normal_standard_alpha_of_mean_difference+inverse_cumulative_normal_standard_alpha_seq)/(1-acceleration*(inverse_cumulative_normal_standard_alpha_of_mean_difference+inverse_cumulative_normal_standard_alpha_seq)))
  #
  # #Accelerated Bootstrap CI
  # Bca1<-stats::quantile(boot_vector, u_adjusted_seq[,1])
  # Bca2<-stats::quantile(boot_vector, u_adjusted_seq[,2])
  # Bca<-cbind(Bca1,Bca2)
  # p.value <- alpha_seq[which.min(0 >= Bca[,1] & 0 <= Bca[,2])]
  #
  # if(p.value==0)
  #   p.value <- 1/( 1 + length(boot_vector))
  #
  # #Adjusted quantiles
  # u_adjusted <- stats::pnorm(inverse_cumulative_normal_standard_alpha_of_mean_difference + (inverse_cumulative_normal_standard_alpha_of_mean_difference+inverse_cumulative_normal_standard_alpha)/(1-acceleration*(inverse_cumulative_normal_standard_alpha_of_mean_difference+inverse_cumulative_normal_standard_alpha)))
  #
  # #Accelerated Bootstrap CI
  # Bca<-stats::quantile(boot_vector, u_adjusted)

  conf.level = 1 - (1 - boot_success/tests_count) * (1 - conf.level)

  Bca <- coxed::bca(boot_vector, conf.level = conf.level)
  p.value <- 0
  if(Bca[1]<0 & Bca[2]>0)
    p.value<- 1


  return(c(Bca, p.value))
}
