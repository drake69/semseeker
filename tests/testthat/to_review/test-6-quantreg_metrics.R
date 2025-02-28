test_that("quantreg_metrics calculates metrics correctly", {

  # Sample data
  predicted_values <- c(1.1, 2.3, 3.0, 4.5)
  expected_values <- c(1.0, 2.0, 3.0, 4.0)
  tau <- 0.5
  res <- data.frame()
  family_test <- "test_family"
  independent_variable <- "test_independent"
  transformation <- "none"
  dependent_variable <- "test_dependent"
  permutation_vector <- c(0.1, 0.2, 0.3, 0.4)

  # Run function without plotting
  result <- quantreg_metrics(predicted_values, expected_values, tau, res, family_test, independent_variable, transformation, dependent_variable, permutation_vector, plot = FALSE)

  # Check the calculated pinball loss
  expect_equal(result$pinball_loss, 0.2)

  # Check the proportion below quantile
  expect_equal(result$below_quantile, 0.75)

  # Check the qq_distance
  expect_equal(result$qq_distance, 0.25)

  # Check the correlation
  expect_equal(result$qq_correlation, cor(expected_values, predicted_values))
})

test_that("quantreg_metrics handles empty data gracefully", {
  # Sample data
  predicted_values <- numeric(0)
  expected_values <- numeric(0)
  tau <- 0.5
  res <- data.frame()
  family_test <- "test_family"
  independent_variable <- "test_independent"
  transformation <- "none"
  dependent_variable <- "test_dependent"
  permutation_vector <- numeric(0)

  # Run function without plotting
  result <- quantreg_metrics(predicted_values, expected_values, tau, res, family_test, independent_variable, transformation, dependent_variable, permutation_vector, plot = FALSE)

  # Expect warnings for empty data
  expect_warning(quantreg_metrics(predicted_values, expected_values, tau, res, family_test, independent_variable, transformation, dependent_variable, permutation_vector, plot = TRUE), "No data available for plotting observed vs. predicted.")
  expect_warning(quantreg_metrics(predicted_values, expected_values, tau, res, family_test, independent_variable, transformation, dependent_variable, permutation_vector, plot = TRUE), "No data available for plotting histogram of betas.")
})

test_that("quantreg_metrics creates plots correctly", {
  # Sample data
  predicted_values <- c(1.1, 2.3, 3.0, 4.5)
  expected_values <- c(1.0, 2.0, 3.0, 4.0)
  tau <- 0.5
  res <- data.frame()
  family_test <- "test_family"
  independent_variable <- "test_independent"
  transformation <- "none"
  dependent_variable <- "test_dependent"
  permutation_vector <- c(0.1, 0.2, 0.3, 0.4)

  # Check that plots are created without error
  expect_error(quantreg_metrics(predicted_values, expected_values, tau, res, family_test, independent_variable, transformation, dependent_variable, permutation_vector, plot = TRUE), NA)
})
