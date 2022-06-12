#' Quantile regression result value, confidence interval and pvalue
#'
#' @param boot_vector vector of boot statistc beta regression
#' @param estimate beta regression
#' @param working_data data to regress
#' @param sig.formula formula for model
#' @param tau quantile to regress at
#' @param independent_variable name of indenpendent variable
#' @param lqm_control controls of lqmm packages
#'
#' @return
quantreg_summary <-function(boot_vector, estimate, working_data, sig.formula, tau, independent_variable, lqm_control){

  alpha_seq <- seq(1e-16, 1-1e-16, 5e-6)
  #BCA method
  #Desired quantiles

  u <- cbind(alpha_seq/2, 1-alpha_seq/2)

  #Compute constants
  z0 <- stats::qnorm(mean(boot_vector <= estimate))
  zu <- stats::qnorm(u)

  # browser()
  n<-nrow(working_data)
  I <- rep(NA, n)
  jack <- rep(NA, n)
  i<-1
  for(i in 1:n){
    #Remove ith working_data point
    xnew <- working_data[-i,]

    #Estimate beta dal modello
    fit.lqm <- suppressMessages(lqmm::lqm( formula = sig.formula, data = xnew, tau = tau, control = lqm_control, fit = TRUE))
    fit.res <- suppressMessages(summary(fit.lqm)$tTable)

    jack[i] <-  fit.res[independent_variable,"Value"]
    I[i] <- (n-1)*(estimate -  jack[i] )
  }
  #Estimate acceleration
  acceleration <- (sum(I^3)/sum(I^2)^1.5)/6

  #Adjusted quantiles
  u_adjusted <- stats::pnorm(z0 + (z0+zu)/(1-acceleration*(z0+zu)))

  #Accelerated Bootstrap CI
  Bca1<-stats::quantile(boot_vector, u_adjusted[,1])
  Bca2<-stats::quantile(boot_vector, u_adjusted[,2])
  selector <- which.min(0 >= Bca1 & 0 <= Bca2)
  pval <- alpha_seq[selector]
  if(pval==0)
    pval <- 1/( 1 + length(boot_vector))


  #Estimate a
  a <- (sum(I^3)/sum(I^2)^1.5)/6
  acc.num <- sum((mean(jack)-jack)^3)
  acc.denom <- sum((mean(jack)-jack)^2)
  accelerate <- acc.num/(6*acc.denom^1.5)
  #Adjusted quantiles
  u_adjusted <- stats::pnorm(z0 + (z0+zu)/(1-a*(z0+zu)))
  #Accelerated Bootstrap CI
  Bca<-stats::quantile(boot_vector, u_adjusted)

  return(c( Bca, pval))
}
