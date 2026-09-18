small_truth <- function(skew = FALSE, nu = NULL) {
  result <- list(
    M = matrix(c(0.2, -0.1, 0.3, 0.4), 2, 2),
    Sigma = matrix(c(1.0, 0.2, 0.2, 0.8), 2, 2),
    Psi = matrix(c(1.0, 0.1, 0.1, 0.9), 2, 2)
  )
  if (skew) result$A <- matrix(c(0.15, -0.05, 0.08, 0.12), 2, 2)
  if (!is.null(nu)) result$nu <- nu
  result
}

expect_positive_definite <- function(X, tolerance = 1e-10) {
  values <- eigen((X + t(X)) / 2, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(is.finite(values)))
  expect_true(all(values > tolerance))
}
