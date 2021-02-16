#' Title
#'
#' @param x vector to compare
#' @param y vector to compare
#'
#' @return true if the order matches otherwise is false

#'
test_match_order <- function(x, y) {
  if (all(x == y)) {
    return(TRUE)
  }
  # print('Perfect match in same order')

  if (!all(x == y) && all(sort(x) == sort(y))) {
    return(FALSE)
  }
  # print('Perfect match in wrong order')

  if (!all(x == y) && !all(sort(x) == sort(y))) {
    return(FALSE)
  }
  # print('No match')
}
