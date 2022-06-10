# bca_pvalue_for_lqm <-
#   function(estimate,
#            boot_vector,
#            model,
#            n_elements)
#   {
#     # browser()
#     # boot.bca <- coxed::bca(na.omit(boot_vector))
#     # p_value <- boot.p_value::boot.p_value(na.omit(boot_vector), type = "bca", theta_null = 0)
#
#     if(stats::shapiro.test(na.omit(boot_vector))$p.value > 0.05)
#     {
#       #boot vector is a normal distribution
#       standard_error <- stats::sd(boot_vector, na.rm = T)
#       mean_boot_vector <- mean(boot_vector, na.rm = T)
#       bias <- estimate - mean_boot_vector
#       bias_threshold <- 0.25 * standard_error
#       p_value <- 2 * stats::pt(-abs((estimate) / standard_error), (n_elements - model$edf))
#       ci_lower_limit <- estimate + stats::qt(0.05 / 2, (n_elements - model$edf)) * standard_error
#       ci_upper_limit <- estimate + stats::qt(1 - 0.05 / 2, (n_elements - model$edf)) * standard_error
#
#       #Bias-corrected p-value
#       if (abs(bias) > bias_threshold) {
#         p_value <- 2 * stats::pt(-abs((2 * estimate - mean_boot_vector) / standard_error), (n_elements - model$edf))
#         ci_lower_limit <-  2 * estimate - mean_boot_vector + stats::qt(0.05 / 2, (n_elements - model$edf)) * standard_error
#         ci_upper_limit <- 2 * estimate - mean_boot_vector + stats::qt(1 - 0.05 / 2, (n_elements - model$edf)) * standard_error
#       }
#       return(c(ci_lower_limit,ci_upper_limit,p_value))
#     }
#     else
#     {
#       # DiCiccio, T. J. and boot_vector. Efron. (1996). Bootstrap Confidence Intervals. Statistical Science. 11(3):
#       # 189–212. https://doi.org/10.1214/ss/1032280214
#       boot.bca <- coxed::bca(stats::na.omit(boot_vector))
#       if(boot.bca[1] <0 & boot.bca[2]>0)
#         p_value <-1
#       else
#         p_value<-0
#       return(c(boot.bca,p_value))
#     }
#
#   }




BCApval<-function(boot_vector, estimate, working_data, sig.formula, tau, independent_variable){

  lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
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
  i<-1
  for(i in 1:n){
    #Remove ith working_data point
    xnew <- working_data[-i,]

    #Estimate beta dal modello
    fit.lqm <- lqmm::lqm( formula = sig.formula, data = xnew, tau = tau, control = lqm_control, fit = TRUE)
    fit.res <- summary(fit.lqm)$tTable

    I[i] <- (n-1)*(estimate -   fit.res[independent_variable,"Value"])
  }
  #Estimate a
  a <- (sum(I^3)/sum(I^2)^1.5)/6

  #Adjusted quantiles
  u_adjusted <- stats::pnorm(z0 + (z0+zu)/(1-a*(z0+zu)))

  #Accelerated Bootstrap CI
  Bca1<-stats::quantile(boot_vector, u_adjusted[,1])
  Bca2<-stats::quantile(boot_vector, u_adjusted[,2])
  Bca<-cbind(Bca1,Bca2)
  selector <- which.min(0 >= Bca[,1] & 0 <= Bca[,2])
  upper <-Bca2[selector]
  lower <- Bca1[selector]
  pval <- alpha_seq[selector]
  if(pval==0)
    pval <- 1/( 1 + length(boot_vector))
  return(c(lower,upper, pval))
}




# semseeker_lqmm<-function(  covariates, independent_variable,burdenValue,family_test,tempDataFrame ){
#
#   # browser()
#   if(is.null(covariates) || length(covariates)==0)
#   {
#     covariates_model <- independent_variable
#   } else
#   {
#     covariates_model <- paste0(paste0(c(independent_variable, covariates),collapse="+", sep=""))
#   }
#   sig.formula <- stats::as.formula(paste0(burdenValue,"~", covariates_model, sep=""))
#   lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
#   quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
#   n_permutations_test <- as.numeric(quantreg_params[3])
#   n_permutations <- as.numeric(quantreg_params[4])
#   tau <- as.numeric(quantreg_params[2])
#
# ####################
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
#     fit.lqm <- lqmm::lqm(model, data = working_data, tau = tau, control = list(verbose = FALSE, loop_tol_ll = 1e-9), fit = TRUE)
#     I[i] <- (n-1)*(estimate -   summary(fit.lqm)[independent_variable,"Value"])
#   }
#   #Estimate a
#   a <- (sum(I^3)/sum(I^2)^1.5)/6
#
#   #Adjusted quantiles
#   u_adjusted <- pnorm(z0 + (z0+zu)/(1-a*(z0+zu)))
#
# ##########################
#
#   model.x <-  lqmm::lqm(sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control)
#   if(n_permutations > n_permutations_test)
#   {
#     model.x.boot <- lqmm::boot(model.x, R = n_permutations_test)
#     beta_full <- summary(model.x.boot)[independent_variable,"Value"]
#     tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
#     colnames(tt) <- colnames(model.x.boot)
#     boot_vector <- stats::na.omit(tt[,independent_variable])
#   }
#
#   #Accelerated Bootstrap CI
#   boot.bca[1]<-quantile(boot_vector, u_adjusted[,1])
#   boot.bca[2]<-quantile(boot_vector, u_adjusted[,2])
#   boot.bca<-cbind(boot.bca[1],boot.bca[2])
#
#   if(boot.bca[1]<0 & boot.bca[2]>0)
#   {
#     n_permutations <- n_permutations_test
#   }
#   else
#   {
#     model.x.boot <- lqmm::boot(model.x, R = n_permutations)
#     beta_full <- summary(model.x.boot)[independent_variable,"Value"]
#     tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
#     colnames(tt) <- colnames(model.x.boot)
#     boot_vector <- stats::na.omit(tt[,independent_variable])
#   }
#
#   #Accelerated Bootstrap CI
#   boot.bca[1]<-quantile(boot_vector, u_adjusted[,1])
#   boot.bca[2]<-quantile(boot_vector, u_adjusted[,2])
#   Bca<-cbind(boot.bca[1],boot.bca[2])
#   pval <- alpha_seq[which.min(0 >= Bca[,1] & 0 <= Bca[,2])]
#   if(pval==0)
#     pval <- 1/( 1 + nrow(tempDataFrame))
#   return(pval)
#
# }




