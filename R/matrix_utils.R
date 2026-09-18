# Canonical numerical utilities used by every model specification.

matrix_vectorize <- function(X) as.vector(X)

frobenius_norm <- function(X) {
  sqrt(sum(as.matrix(X)^2))
}

relative_frobenius_error <- function(estimate, truth) {
  denominator <- frobenius_norm(truth)
  if (denominator <= .Machine$double.eps) {
    return(frobenius_norm(estimate - truth))
  }
  frobenius_norm(estimate - truth) / denominator
}

is_symmetric_matrix <- function(X, tolerance = 1e-8) {
  is.matrix(X) && nrow(X) == ncol(X) &&
    isTRUE(all.equal(X, t(X), tolerance = tolerance))
}

matrix_is_posdef <- function(X, precision = 1e-10) {
  X <- as.matrix(X)
  if (!is_symmetric_matrix(X, tolerance = sqrt(precision)) ||
      any(!is.finite(X))) return(FALSE)
  values <- eigen(matrix_symmetrize(X), symmetric = TRUE,
                  only.values = TRUE)$values
  all(values > precision)
}

matrix_assert_posdef <- function(X, name, precision = 1e-10) {
  if (!matrix_is_posdef(X, precision)) {
    stop(sprintf("%s must be symmetric positive definite.", name), call. = FALSE)
  }
  invisible(TRUE)
}

matrix_symmetrize <- function(X) {
  X <- as.matrix(X)
  (X + t(X)) / 2
}

matrix_make_posdef <- function(X, epsilon = 1e-8, eig_floor = NULL) {
  if (!is.null(eig_floor)) epsilon <- eig_floor
  X <- as.matrix(X)
  if (nrow(X) != ncol(X)) {
    stop("The matrix to be regularized must be square.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("The matrix to be regularized contains non-finite values.", call. = FALSE)
  }
  X <- matrix_symmetrize(X)
  decomposition <- eigen(X, symmetric = TRUE)
  values <- pmax(decomposition$values, epsilon)
  result <- decomposition$vectors %*%
    diag(values, nrow = length(values)) %*%
    t(decomposition$vectors)
  matrix_symmetrize(result)
}

matrix_chol_inverse <- function(X, epsilon = 1e-8) {
  chol2inv(chol(matrix_make_posdef(X, epsilon)))
}

matrix_normalize_Psi <- function(Psi, epsilon = 1e-8) {
  Psi <- matrix_make_posdef(Psi, epsilon)
  log_determinant <- as.numeric(determinant(Psi, logarithm = TRUE)$modulus)
  Psi / exp(log_determinant / ncol(Psi))
}

ecm_relative_change <- function(current, previous) {
  if (!is.finite(current) || !is.finite(previous) || current == 0) return(Inf)
  abs(1 - previous / current)
}

model_parameter_count <- function(p, q, skew = FALSE, extra = 0L) {
  p * q + as.integer(skew) * p * q +
    p * (p + 1) / 2 + q * (q + 1) / 2 - 1 + extra
}

model_bic <- function(loglik, parameter_count, n) {
  -2 * loglik + parameter_count * log(n)
}
