substitute_infinite <- function(x) {
  # browser()
  x <- as.matrix(x)
  inf_vals <- !is.finite(x)
  # get row and column of infinite values
  inf_values <- which(inf_vals, arr.ind = TRUE)
  # count the number of infinite values
  n_inf_values <- nrow(inf_values)
  # get max abs value
  max_abs_value <- max(abs(x[is.finite(x)]), na.rm = TRUE)
  # if there are infinite values
  if (n_inf_values > 0) {
    # for each infinite value
    for (i in 1:n_inf_values) {
      # get the row and column of the infinite value
      row <- inf_values[i, 1]
      col <- inf_values[i, 2]
      # get the value of the infinite value
      value <- x[row, col]
      # if the value is infinite
      if (!is.finite(value)) {
        # replace the infinite value with max_abs and correct sign
        x[row, col] <- sign(value) * max_abs_value
      }
    }
  }
  # table(is.infinite(x))
  as.data.frame(x)
}
