#' Quantile regression result value, confidence interval and pvalue
#'
#' @param boot_vector vector of boot statistc beta regression
#' @param estimate beta regression
#' @param conf.level confidence intervals alpha level
#' @param boot_success number of success respecting the null hypothesis
#' @param tests_count how many tests were done
#'
#' @return ci and pvalue with BCA method
#' @importFrom doRNG %dorng%
quantreg_summary <-function(boot_vector, estimate, conf.level , boot_success = 0, tests_count=1 ){

  conf.level = 1 - (1 - boot_success/tests_count) * (1 - conf.level)
  Bca <- coxed::bca(boot_vector, conf.level = conf.level)
  p.value <- mean(abs(boot_vector) >= abs(estimate))

  return(c(Bca, p.value))
}
