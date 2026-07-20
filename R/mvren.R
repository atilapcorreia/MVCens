# Matrix-Variate Row Exponential-Normal (MVREN)
#
# Self-contained implementation of the model
#
#   X = M + W A + V,
#
# where W = diag(W_1, ..., W_p), W_i are mutually independent with
# W_i ~ Exp(lambda_i), and V ~ N_{p x q}(0, Sigma, Psi).
#
# Main functions:
#   dmvren()              density
#   loglik_mvren()        observed-data log-likelihood
#   rmvren()              random generation
#   mvren_ecm()           ECM estimation
#   mvren_mean()          theoretical mean
#   mvren_covariances()   row/column covariance summaries
#   mvren_monte_carlo()   simulation audit
#
#   Required packages:
#   mvtnorm, tmvtnorm

# ---------------------------------------------------------------------------
# Internal numerical utilities
# ---------------------------------------------------------------------------

#' Vectorize a matrix
#'
#' Converts a matrix or matrix-like object into a vector using the default
#' column-major ordering employed by R.
#'
#' @param B A matrix or matrix-like object to be vectorized.
#'
#' @return A vector containing the elements of `B` arranged column by column.
#'
#' @keywords internal
#' @noRd
mvren_vec <- function(B) {
  as.vector(B)
}

#' Compute the Frobenius norm
#'
#' Computes the Frobenius norm of a matrix or matrix-like object by taking
#' the square root of the sum of its squared elements.
#'
#' @param X A matrix or matrix-like object.
#'
#' @return A nonnegative numeric scalar representing the Frobenius norm of `X`.
#'
#' @keywords internal
#' @noRd
mvren_frobenius <- function(X) {
  sqrt(sum(as.matrix(X)^2))
}

#' Compute the relative Frobenius error
#'
#' Computes the Frobenius norm of the difference between an estimated object
#' and its true value, normalized by the Frobenius norm of the true value.
#' When the norm of the true value is numerically zero, the unnormalized
#' Frobenius error is returned.
#'
#' @param estimate A matrix or matrix-like object containing the estimated
#'   values.
#' @param truth A matrix or matrix-like object containing the reference or
#'   true values.
#'
#' @return A nonnegative numeric scalar representing the relative Frobenius
#'   error, or the absolute Frobenius error when `truth` has a numerically
#'   zero norm.
#'
#' @keywords internal
#' @noRd
mvren_relative_frobenius <- function(estimate, truth) {
  denom <- mvren_frobenius(truth)
  if (denom <= .Machine$double.eps) {
    return(mvren_frobenius(estimate - truth))
  }
  mvren_frobenius(estimate - truth) / denom
}
#' Check matrix symmetry
#'
#' Determines whether an object is a square matrix that is numerically
#' symmetric within a specified tolerance. The comparison is performed with
#' [all.equal()], so matrix attributes must also be compatible.
#'
#' @param X An object to be tested for symmetry.
#' @param tol A nonnegative numeric scalar specifying the tolerance used when
#'   comparing `X` with its transpose. The default is `1e-8`.
#'
#' @return A logical scalar indicating whether `X` is a square matrix and is
#'   symmetric, with compatible attributes, within the specified tolerance.
#'
#' @keywords internal
#' @noRd
mvren_is_symmetric <- function(X, tol = 1e-8) {
  is.matrix(X) &&
    nrow(X) == ncol(X) &&
    isTRUE(all.equal(X, t(X), tolerance = tol))
}

#' Check positive definiteness
#'
#' Determines whether an object is a numerically symmetric positive-definite
#' matrix by examining its eigenvalues.
#'
#' @param X A matrix or matrix-like object to be tested.
#' @param tol A positive numeric scalar specifying the minimum admissible
#'   eigenvalue. The default is `1e-10`.
#'
#' @return A logical scalar indicating whether `X` is symmetric within the
#'   numerical tolerance and all of its eigenvalues are finite and greater
#'   than `tol`.
#'
#' @keywords internal
#' @noRd
mvren_is_posdef <- function(X, tol = 1e-10) {
  X <- as.matrix(X)
  if (!mvren_is_symmetric(X, tol = sqrt(tol))) {
    return(FALSE)
  }
  values <- eigen((X + t(X)) / 2, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(values)) && all(values > tol)
}

#' Assert positive definiteness
#'
#' Verifies that an object is a symmetric positive-definite matrix within a
#' specified numerical tolerance. An error is raised when the condition is
#' not satisfied.
#'
#' @param X A matrix or matrix-like object to be checked.
#' @param name A character string identifying `X` in the error message.
#' @param tol A positive numeric scalar specifying the minimum admissible
#'   eigenvalue. The default is `1e-10`.
#'
#' @return Invisibly returns `TRUE` when `X` is symmetric positive definite.
#'   Otherwise, the function terminates with an error.
#'
#' @keywords internal
#' @noRd
mvren_assert_posdef <- function(X, name, tol = 1e-10) {
  if (!mvren_is_posdef(X, tol = tol)) {
    stop(sprintf("%s must be symmetric positive definite.", name), call. = FALSE)
  }
  invisible(TRUE)
}
#' Regularize a matrix using a positive eigenvalue floor
#'
#' Symmetrizes a square finite matrix and replaces every eigenvalue smaller
#' than `eig_floor` with `eig_floor`.
#'
#' @param S A square numeric matrix or matrix-like object to be regularized.
#' @param eig_floor A numeric scalar specifying the minimum permitted
#'   eigenvalue. The default is `1e-8`. A strictly positive value is required
#'   for the result to be positive definite in exact arithmetic.
#'
#' @return A symmetric matrix whose eigenvalues are bounded below by
#'   `eig_floor`. The result is positive definite when `eig_floor > 0`, subject
#'   to numerical rounding; when `eig_floor = 0`, it may be only positive
#'   semidefinite.
#'
#' @keywords internal
#' @noRd
mvren_make_posdef <- function(S, eig_floor = 1e-8) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S)) {
    stop("The matrix to be regularized must be square.", call. = FALSE)
  }
  if (any(!is.finite(S))) {
    stop("The matrix to be regularized contains non-finite values.", call. = FALSE)
  }

  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  values <- pmax(ev$values, eig_floor)
  S_pd <- ev$vectors %*% diag(values, nrow = length(values)) %*% t(ev$vectors)
  (S_pd + t(S_pd)) / 2
}
#' Prepare a numerical precision matrix
#'
#' Symmetrizes a matrix and determines whether its smallest eigenvalue is at
#' least `eig_floor`. When the smallest eigenvalue falls below this threshold,
#' a diagonal correction is handled according to the selected policy.
#'
#' @param Q A square numeric matrix or matrix-like object to be checked and
#'   optionally regularized.
#' @param policy A character string specifying how to handle a matrix whose
#'   smallest eigenvalue is below `eig_floor`. Available options are `"warn"`,
#'   which regularizes the matrix and issues a warning; `"strict"`, which
#'   terminates with an error; and `"regularize"`, which regularizes the matrix
#'   silently.
#' @param eig_floor A numeric scalar specifying the minimum admissible
#'   eigenvalue. The default is `1e-8`.
#'
#' @return A list containing:
#' \describe{
#'   \item{Q}{The symmetrized and, when necessary, regularized matrix.}
#'   \item{min_eigenvalue}{The smallest eigenvalue before regularization.}
#'   \item{regularization}{The diagonal correction applied to the matrix.}
#'   \item{was_regularized}{A logical scalar indicating whether regularization
#'   was applied.}
#' }
#'
#' @keywords internal
#' @noRd
mvren_prepare_Q <- function(Q, policy = c("warn", "strict", "regularize"),
                            eig_floor = 1e-8) {
  policy <- match.arg(policy)
  Q <- (as.matrix(Q) + t(as.matrix(Q))) / 2
  values <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(values)
  correction <- max(eig_floor - minimum, 0)
  if (correction > 0 && policy == "strict") {
    stop(sprintf("Q is not numerically positive definite (min eigenvalue=%g).",
                 minimum), call. = FALSE)
  }
  if (correction > 0 && policy == "warn") {
    warning(sprintf("Q regularized by a diagonal correction of %g.", correction),
            call. = FALSE)
  }
  if (correction > 0) Q <- Q + correction * diag(nrow(Q))
  list(Q = Q, min_eigenvalue = minimum, regularization = correction,
       was_regularized = correction > 0)
}
#' Regularize a matrix to be positive semidefinite
#'
#' Symmetrizes a square finite matrix and replaces eigenvalues smaller than
#' `eig_floor` with that threshold. This low-level helper assumes that `S` is
#' square and finite and that `eig_floor` is finite and nonnegative.
#'
#' @param S A square numeric matrix or matrix-like object to be regularized.
#' @param eig_floor A nonnegative numeric scalar specifying the minimum
#'   permitted eigenvalue. The default is `0`.
#'
#' @return A symmetric positive-semidefinite matrix whose eigenvalues are
#'   bounded below by `eig_floor`. If `eig_floor > 0`, the returned matrix is
#'   positive definite in exact arithmetic, subject to numerical rounding.
#'
#' @keywords internal
#' @noRd
mvren_make_psd <- function(S, eig_floor = 0) {
  S <- as.matrix(S)
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  values <- pmax(ev$values, eig_floor)
  S_psd <- ev$vectors %*% diag(values, nrow = length(values)) %*% t(ev$vectors)
  (S_psd + t(S_psd)) / 2
}

#' Compute a lower-triangular Cholesky factor
#'
#' Computes the lower-triangular Cholesky factor of a symmetric
#' positive-definite matrix. The input matrix is validated before the
#' factorization is performed.
#'
#' @param S A symmetric positive-definite matrix or matrix-like object.
#' @param name A character string identifying `S` in any error message.
#'   The default is `"matrix"`.
#'
#' @return A lower-triangular matrix `L` satisfying
#'   `S = L %*% t(L)`.
#'
#' @keywords internal
#' @noRd
mvren_chol_lower <- function(S, name = "matrix") {
  S <- as.matrix(S)
  mvren_assert_posdef(S, name)
  t(chol(S))
}

#' Solve a linear system using a Cholesky factorization
#'
#' Solves a linear system with a symmetric positive-definite coefficient
#' matrix using its lower-triangular Cholesky factor. The solution is obtained
#' through successive forward and backward substitutions.
#'
#' @param S A symmetric positive-definite coefficient matrix or matrix-like
#'   object.
#' @param B A numeric vector or matrix representing the right-hand side of the
#'   linear system.
#' @param name A character string identifying `S` in any error message.
#'   The default is `"matrix"`.
#'
#' @return A numeric vector or matrix containing the solution to
#'   `S %*% X = B`.
#'
#' @keywords internal
#' @noRd
mvren_chol_solve <- function(S, B, name = "matrix") {
  L <- mvren_chol_lower(S, name)
  backsolve(t(L), forwardsolve(L, B))
}

#' Compute the logarithm of a matrix determinant
#'
#' Computes the natural logarithm of the determinant of a matrix using a
#' numerically stable logarithmic determinant calculation. An error is raised
#' when the determinant is nonpositive.
#'
#' @param S A square numeric matrix or matrix-like object.
#' @param name A character string identifying `S` in any error message.
#'   The default is `"matrix"`.
#'
#' @return A numeric scalar containing the natural logarithm of the
#'   determinant of `S`.
#'
#' @keywords internal
#' @noRd
mvren_logdet <- function(S, name = "matrix") {
  det_obj <- determinant(S, logarithm = TRUE)
  if (det_obj$sign <= 0) {
    stop(sprintf("%s must have a positive determinant.", name), call. = FALSE)
  }
  as.numeric(det_obj$modulus)
}
#' Validate the rate parameter vector
#'
#' Coerces the supplied object to a numeric vector and verifies that it has the
#' required length and contains only finite, strictly positive values.
#'
#' @param lambda An object coercible to a numeric vector containing the rate
#'   parameters.
#' @param p A positive integer specifying the required length of `lambda`.
#'
#' @return The validated rate parameter vector as a numeric vector.
#'
#' @keywords internal
#' @noRd
mvren_validate_lambda <- function(lambda, p) {
  lambda <- as.numeric(lambda)
  if (length(lambda) != p) {
    stop(sprintf("lambda must be a numeric vector of length %d.", p), call. = FALSE)
  }
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("Every component of lambda must be finite and strictly positive.", call. = FALSE)
  }
  lambda
}

#' Validate MVREN model parameters
#'
#' Validates the dimensions, numeric types, rate parameters, and, optionally,
#' the positive definiteness of the covariance matrices used in the MVREN
#' model. The inputs are converted to matrix or numeric-vector form before
#' validation.
#'
#' @param X An optional numeric matrix representing an observation. When
#'   supplied, it must have the same dimensions as `M` and `A`.
#' @param M A numeric `p` by `q` location matrix.
#' @param A A numeric `p` by `q` latent-effect matrix.
#' @param Sigma A numeric `p` by `p` row covariance matrix.
#' @param Psi A numeric `q` by `q` column covariance matrix.
#' @param lambda An optional numeric vector of length `p` containing finite,
#'   strictly positive rate parameters. When `NULL`, a vector of ones is used.
#' @param check_posdef A logical scalar indicating whether `Sigma` and `Psi`
#'   should be checked for positive definiteness. The default is `TRUE`.
#'
#' @return A list containing the validated matrices `M`, `A`, `Sigma`, and
#'   `Psi`; the validated rate vector `lambda`; and the matrix dimensions `p`
#'   and `q`.
#'
#' @keywords internal
#' @noRd
mvren_validate_parameters <- function(X = NULL,
                                      M,
                                      A,
                                      Sigma,
                                      Psi,
                                      lambda = NULL,
                                      check_posdef = TRUE) {
  M <- as.matrix(M)
  A <- as.matrix(A)
  Sigma <- as.matrix(Sigma)
  Psi <- as.matrix(Psi)

  if (!is.numeric(M) || !is.numeric(A) || !is.numeric(Sigma) || !is.numeric(Psi)) {
    stop("M, A, Sigma, and Psi must be numeric matrices.", call. = FALSE)
  }

  p <- nrow(M)
  q <- ncol(M)

  if (!all(dim(A) == c(p, q))) {
    stop("A must have the same dimensions as M.", call. = FALSE)
  }
  if (!all(dim(Sigma) == c(p, p))) {
    stop("Sigma must be p x p, where p = nrow(M).", call. = FALSE)
  }
  if (!all(dim(Psi) == c(q, q))) {
    stop("Psi must be q x q, where q = ncol(M).", call. = FALSE)
  }

  if (!is.null(X)) {
    X <- as.matrix(X)
    if (!all(dim(X) == c(p, q))) {
      stop("X must have the same dimensions as M and A.", call. = FALSE)
    }
  }

  if (is.null(lambda)) {
    lambda <- rep(1, p)
  }
  lambda <- mvren_validate_lambda(lambda, p)

  if (check_posdef) {
    mvren_assert_posdef(Sigma, "Sigma")
    mvren_assert_posdef(Psi, "Psi")
  }

  list(
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda,
    p = p,
    q = q
  )
}

#' Validate an MVREN sample array
#'
#' Verifies that the supplied sample is a three-dimensional numeric array with
#' dimensions `p` by `q` by `n`, strictly positive dimensions, and only finite
#' values.
#'
#' @param X An array-like object representing a sample of matrix-valued
#'   observations.
#'
#' @return The validated sample as a three-dimensional numeric array.
#'
#' @keywords internal
#' @noRd
mvren_validate_sample_array <- function(X) {
  X <- as.array(X)
  dims <- dim(X)
  if (length(dims) != 3L) {
    stop("X must be a 3-dimensional array with dimensions p x q x n.", call. = FALSE)
  }
  if (any(dims <= 0L)) {
    stop("All dimensions of X must be positive.", call. = FALSE)
  }
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must contain only finite numeric values.", call. = FALSE)
  }
  X
}

#' Require an installed package namespace
#'
#' Verifies that a package namespace is available without attaching the
#' package to the search path. An informative error is raised when the package
#' is not installed.
#'
#' @param package A character string specifying the name of the required
#'   package.
#'
#' @return Invisibly returns `TRUE` when the package namespace is available.
#'   Otherwise, the function terminates with an error containing installation
#'   instructions.
#'
#' @keywords internal
#' @noRd
mvren_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required. Install it with install.packages('%s').", package, package),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
#' Compute a log multivariate normal cumulative probability
#'
#' Computes the logarithm of the lower-orthant multivariate normal probability
#' that every component is less than or equal to its corresponding upper
#' integration limit. The probability is evaluated using \pkg{mvtnorm}.
#'
#' @param upper A numeric vector specifying the upper integration limits.
#' @param mean A numeric vector specifying the mean of the multivariate normal
#'   distribution.
#' @param sigma A symmetric positive-definite covariance matrix.
#'
#' @return A numeric scalar containing the logarithm of the multivariate normal
#'   cumulative probability. Returns `-Inf` when the computed probability is
#'   non-finite or nonpositive.
#'
#' @keywords internal
#' @noRd
mvren_log_pmvnorm_upper <- function(upper, mean, sigma) {
  mvren_require_namespace("mvtnorm")
  prob <- suppressWarnings(
    as.numeric(mvtnorm::pmvnorm(upper = upper, mean = mean, sigma = sigma))
  )
  if (!is.finite(prob) || prob <= 0) {
    return(-Inf)
  }
  log(prob)
}
#' Compute positive-orthant truncated-normal moments
#'
#' Computes the conditional first and second moments of a multivariate normal
#' distribution truncated to the positive orthant. Closed-form componentwise
#' formulas are used in the univariate case and when the covariance matrix is
#' diagonal. For a general covariance matrix, the moments are evaluated using
#' \pkg{tmvtnorm}.
#'
#' @param mean A numeric vector specifying the mean of the untruncated normal
#'   distribution.
#' @param sigma A covariance matrix with dimensions compatible with `mean`.
#' @param eig_floor A positive numeric scalar used as a lower bound for
#'   marginal variances in the univariate and diagonal calculations. When a
#'   variance is below this threshold, the returned moments correspond to the
#'   stabilized variance rather than the original value. The default is
#'   `1e-8`.
#' @param max_retries A positive integer specifying the maximum number of
#'   repeated attempts to obtain valid multivariate truncated-normal moments.
#'   The default is `4L`.
#'
#' @return A list containing:
#' \describe{
#'   \item{mean}{A numeric vector containing the conditional mean under
#'   positive-orthant truncation.}
#'   \item{variance}{The conditional covariance matrix under
#'   positive-orthant truncation.}
#'   \item{second_moment}{The conditional second-moment matrix, equal to the
#'   conditional covariance matrix plus the outer product of the conditional
#'   mean.}
#' }
#'
#' @keywords internal
#' @noRd
mvren_truncated_moments <- function(mean,
                                    sigma,
                                    eig_floor = 1e-8,
                                    max_retries = 4L) {
  mean <- as.numeric(mean)
  p <- length(mean)
  sigma <- as.matrix(sigma)

  if (!all(dim(sigma) == c(p, p))) {
    stop("The truncated-normal covariance has incompatible dimensions.", call. = FALSE)
  }

  univariate_moments <- function(location, variance) {
    variance <- max(as.numeric(variance), eig_floor)
    sd_value <- sqrt(variance)
    alpha <- -location / sd_value
    log_tail <- stats::pnorm(alpha, lower.tail = FALSE, log.p = TRUE)
    ratio <- exp(stats::dnorm(alpha, log = TRUE) - log_tail)
    truncated_mean <- location + sd_value * ratio
    truncated_variance <- variance * (1 + alpha * ratio - ratio^2)
    truncated_variance <- max(truncated_variance, .Machine$double.eps)
    c(mean = truncated_mean, variance = truncated_variance)
  }

  # Exact univariate moments are more stable than a generic multivariate call.
  if (p == 1L) {
    moments <- univariate_moments(mean[1], sigma[1, 1])
    return(list(
      mean = unname(moments["mean"]),
      variance = matrix(unname(moments["variance"]), 1, 1),
      second_moment = matrix(unname(moments["variance"] + moments["mean"]^2), 1, 1)
    ))
  }

  off_diagonal <- sigma
  diag(off_diagonal) <- 0
  diagonal_tolerance <- 100 * .Machine$double.eps * max(1, max(abs(sigma)))
  if (max(abs(off_diagonal)) <= diagonal_tolerance) {
    # With diagonal covariance, truncation by the positive orthant preserves
    # independence.  Componentwise moments are exact and avoid multivariate
    # tail underflow and spurious off-diagonal covariances from quadrature.
    moments <- vapply(seq_len(p), function(index) {
      univariate_moments(mean[index], sigma[index, index])
    }, numeric(2))
    truncated_mean <- unname(moments["mean", ])
    truncated_variance <- diag(unname(moments["variance", ]), p, p)
    return(list(
      mean = truncated_mean,
      variance = truncated_variance,
      second_moment = truncated_variance + tcrossprod(truncated_mean)
    ))
  }

  mvren_require_namespace("tmvtnorm")
  last_error <- NULL

  for (attempt in seq_len(max_retries)) {
    result <- tryCatch(
      tmvtnorm::mtmvnorm(
        mean = mean,
        sigma = sigma,
        lower = rep(0, p),
        upper = rep(Inf, p)
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      last_error <- conditionMessage(result)
      next
    }

    truncated_mean <- as.numeric(result$tmean)
    truncated_variance <- as.matrix(result$tvar)

    valid_mean <- length(truncated_mean) == p && all(is.finite(truncated_mean))
    valid_variance <- all(dim(truncated_variance) == c(p, p)) &&
      all(is.finite(truncated_variance))

    if (!valid_mean || !valid_variance) {
      last_error <- "tmvtnorm returned invalid conditional moments."
      next
    }

    truncated_variance <- mvren_make_psd(truncated_variance, eig_floor = 0)
    second_moment <- truncated_variance + tcrossprod(truncated_mean)

    return(list(
      mean = truncated_mean,
      variance = truncated_variance,
      second_moment = second_moment
    ))
  }

  stop(
    paste0(
      "Unable to compute valid positive-orthant truncated-normal moments after ",
      max_retries,
      " attempts. Last error: ",
      if (is.null(last_error)) "unknown numerical failure" else last_error
    ),
    call. = FALSE
  )
}

#' Apply an identifiability constraint to A and lambda
#'
#' Applies the selected identifiability convention to the latent-effect matrix
#' and rate vector. Under the unit-row-norm constraint, each row of `A` is
#' normalized to have Euclidean norm one and the corresponding component of
#' `lambda` is rescaled to preserve the distribution of the latent
#' contribution.
#'
#' @param A A numeric matrix whose rows contain the latent-effect directions.
#' @param lambda A numeric vector containing one finite, strictly positive rate
#'   parameter for each row of `A`.
#' @param mode A character string specifying the identifiability convention.
#'   Available options are `"fixed"`, `"unit_A_rows"`, and `"unconstrained"`.
#'   Only `"unit_A_rows"` modifies `A` and `lambda`.
#' @param row_norm_floor A positive numeric scalar specifying the minimum
#'   admissible row norm of `A` when `mode = "unit_A_rows"`. The default is
#'   `1e-10`.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{The original or row-normalized latent-effect matrix.}
#'   \item{lambda}{The original or correspondingly rescaled rate vector.}
#' }
#'
#' @keywords internal
#' @noRd
mvren_identify_A_lambda <- function(A,
                                    lambda,
                                    mode = c("fixed", "unit_A_rows", "unconstrained"),
                                    row_norm_floor = 1e-10) {
  mode <- match.arg(mode)
  A <- as.matrix(A)
  lambda <- mvren_validate_lambda(lambda, nrow(A))

  if (mode != "unit_A_rows") {
    return(list(A = A, lambda = lambda))
  }

  row_norms <- sqrt(rowSums(A^2))
  if (any(!is.finite(row_norms)) || any(row_norms <= row_norm_floor)) {
    stop(
      paste0(
        "A row is numerically zero, so lambda is not identifiable under the ",
        "unit-row-norm constraint. Use lambda_mode = 'fixed' or provide better starting values."
      ),
      call. = FALSE
    )
  }

  # The transformation (A_i, lambda_i) -> (A_i/c_i, lambda_i/c_i),
  # with c_i = ||A_i||, preserves the distribution of W_i A_i while imposing
  # ||A_i||_2 = 1.
  list(
    A = A / row_norms,
    lambda = lambda / row_norms
  )
}

#' Construct initial values for the MVREN model
#'
#' Constructs starting values for the MVREN model parameters. When the
#' latent-effect matrix is not supplied, its row directions are initialized
#' from elementwise empirical skewness coefficients and scaled according to
#' `skew_scale`. Default covariance matrices are identity matrices, and the
#' default location matrix is chosen so that the initial model mean matches
#' the sample mean.
#'
#' @param X A three-dimensional numeric array with dimensions `p` by `q` by
#'   `n`, containing the matrix-valued observations.
#' @param M_init An optional numeric `p` by `q` matrix providing the initial
#'   value of the location matrix. When `NULL`, it is constructed so that the
#'   initial model mean equals the sample mean.
#' @param A_init An optional numeric `p` by `q` matrix providing the initial
#'   value of the latent-effect matrix. When `NULL`, it is initialized using
#'   elementwise empirical skewness coefficients.
#' @param Sigma_init An optional symmetric positive-definite `p` by `p` matrix
#'   providing the initial row covariance matrix. When `NULL`, the identity
#'   matrix is used.
#' @param Psi_init An optional symmetric positive-definite `q` by `q` matrix
#'   providing the initial column covariance matrix. When `NULL`, the identity
#'   matrix is used.
#' @param lambda_init An optional numeric vector of length `p` containing
#'   finite, strictly positive initial rate parameters. When `NULL`, a vector
#'   of ones is used.
#' @param lambda_mode A character string specifying the identifiability
#'   convention for `A` and `lambda`. Available options are `"fixed"`,
#'   `"unit_A_rows"`, and `"unconstrained"`. Under `"unit_A_rows"`, the rows
#'   of `A` are normalized and `lambda` is correspondingly rescaled.
#' @param skew_scale A positive numeric scalar controlling the Euclidean norm
#'   assigned to each automatically initialized row of `A`. The default is
#'   `0.5`.
#'
#' @return A list containing the initial parameter values:
#' \describe{
#'   \item{M}{The initial `p` by `q` location matrix.}
#'   \item{A}{The initial `p` by `q` latent-effect matrix.}
#'   \item{Sigma}{The initial `p` by `p` row covariance matrix.}
#'   \item{Psi}{The initial `q` by `q` column covariance matrix.}
#'   \item{lambda}{The initial rate parameter vector of length `p`.}
#' }
#'
#' @keywords internal
#' @noRd
mvren_initial_values <- function(X,
                                 M_init = NULL,
                                 A_init = NULL,
                                 Sigma_init = NULL,
                                 Psi_init = NULL,
                                 lambda_init = NULL,
                                 lambda_mode = c("fixed", "unit_A_rows", "unconstrained"),
                                 skew_scale = 0.5) {
  X <- mvren_validate_sample_array(X)
  lambda_mode <- match.arg(lambda_mode)

  p <- dim(X)[1]
  q <- dim(X)[2]
  n <- dim(X)[3]

  lambda <- if (is.null(lambda_init)) rep(1, p) else mvren_validate_lambda(lambda_init, p)
  Xbar <- apply(X, c(1, 2), mean)

  if (is.null(A_init)) {
    centered <- X
    for (i in seq_len(n)) {
      centered[, , i] <- X[, , i] - Xbar
    }

    sds <- apply(X, c(1, 2), stats::sd)
    sds[!is.finite(sds) | sds <= .Machine$double.eps] <- 1
    skew <- apply(centered, c(1, 2), function(z) mean(z^3)) / (sds^3)
    skew[!is.finite(skew)] <- 0
    A <- matrix(0, p, q)
    for (row in seq_len(p)) {
      direction <- skew[row, ]
      direction_norm <- sqrt(sum(direction^2))
      if (!is.finite(direction_norm) ||
          direction_norm <= .Machine$double.eps^0.25) {
        direction <- ((seq_len(q) + row) %% 2) * 2 - 1
        direction_norm <- sqrt(sum(direction^2))
      }
      A[row, ] <- skew_scale * direction / direction_norm
    }
  } else {
    A <- as.matrix(A_init)
  }

  if (!all(dim(A) == c(p, q))) {
    stop("A_init must have dimensions p x q.", call. = FALSE)
  }

  identified <- mvren_identify_A_lambda(A, lambda, mode = lambda_mode)
  A <- identified$A
  lambda <- identified$lambda

  M <- if (is.null(M_init)) {
    Xbar - diag(1 / lambda, p, p) %*% A
  } else {
    as.matrix(M_init)
  }

  Sigma <- if (is.null(Sigma_init)) diag(p) else as.matrix(Sigma_init)
  Psi <- if (is.null(Psi_init)) diag(q) else as.matrix(Psi_init)

  mvren_validate_parameters(
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda
  )

  list(M = M, A = A, Sigma = Sigma, Psi = Psi, lambda = lambda)
}


# ---------------------------------------------------------------------------
# Theoretical moments
# ---------------------------------------------------------------------------

#' Compute the mean matrix of an MVREN distribution
#'
#' Computes the mean matrix of a matrix-variate row exponential-normal
#' distribution from its location matrix, latent-effect matrix, and
#' row-specific exponential rate parameters.
#'
#' @param M A numeric `p` by `q` location matrix.
#' @param A A numeric `p` by `q` latent-effect matrix.
#' @param lambda An optional numeric vector of length `p` containing finite,
#'   strictly positive exponential rate parameters. When `NULL`, a vector of
#'   ones is used.
#'
#' @return A numeric `p` by `q` matrix containing the mean of the MVREN
#'   distribution.
#'
#' @export
mvren_mean <- function(M, A, lambda = NULL) {
  M <- as.matrix(M)
  A <- as.matrix(A)
  if (!all(dim(M) == dim(A))) {
    stop("M and A must have the same dimensions.", call. = FALSE)
  }
  p <- nrow(M)
  if (is.null(lambda)) {
    lambda <- rep(1, p)
  }
  lambda <- mvren_validate_lambda(lambda, p)
  M + diag(1 / lambda, p, p) %*% A
}
#' Compute row and column covariance summaries of an MVREN distribution
#'
#' Computes the row and column covariance summaries associated with a
#' matrix-variate row exponential-normal distribution from its latent-effect
#' matrix, row and column covariance matrices, and row-specific exponential
#' rate parameters.
#'
#' @param M A numeric `p` by `q` location matrix. The covariance summaries do
#'   not depend on `M`; it is used to establish and validate the model
#'   dimensions.
#' @param A A numeric `p` by `q` latent-effect matrix.
#' @param Sigma A symmetric positive-definite `p` by `p` row covariance
#'   matrix.
#' @param Psi A symmetric positive-definite `q` by `q` column covariance
#'   matrix.
#' @param lambda An optional numeric vector of length `p` containing finite,
#'   strictly positive exponential rate parameters. When `NULL`, a vector of
#'   ones is used.
#'
#' @return A list containing:
#' \describe{
#'   \item{row}{A \eqn{p \times p} matrix equal to
#'   \eqn{\mathbb{E}\left[
#'   (\mathbf{X}-\mathbb{E}(\mathbf{X}))
#'   (\mathbf{X}-\mathbb{E}(\mathbf{X}))^{\top}
#'   \right]}.}
#'   \item{column}{A \eqn{q \times q} matrix equal to
#'   \eqn{\mathbb{E}\left[
#'   (\mathbf{X}-\mathbb{E}(\mathbf{X}))^{\top}
#'   (\mathbf{X}-\mathbb{E}(\mathbf{X}))
#'   \right]}.}
#' }
#'
#' @export
mvren_covariances <- function(M, A, Sigma, Psi, lambda = NULL) {
  pars <- mvren_validate_parameters(
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda
  )

  latent_variances <- 1 / pars$lambda^2
  AAt <- pars$A %*% t(pars$A)

  list(
    row = diag(latent_variances * diag(AAt), pars$p, pars$p) +
      sum(diag(pars$Psi)) * pars$Sigma,
    column = t(pars$A) %*% diag(latent_variances, pars$p, pars$p) %*% pars$A +
      sum(diag(pars$Sigma)) * pars$Psi
  )
}


#' Density of the matrix-variate row exponential-normal distribution
#'
#' Computes the density or log-density of the matrix-variate row
#' exponential-normal (MVREN) distribution at a matrix-valued observation.
#' The MVREN model extends the matrix-variate normal distribution by
#' introducing independent, nonnegative row-specific latent effects, allowing
#' different rows of the observation matrix to exhibit distinct directions
#' and degrees of asymmetry.
#'
#' @param X A numeric `p` by `q` observation matrix.
#' @param M A numeric `p` by `q` location matrix.
#' @param A A numeric `p` by `q` row-specific latent-effect matrix.
#' @param Sigma A symmetric positive-definite `p` by `p` row covariance
#'   matrix.
#' @param Psi A symmetric positive-definite `q` by `q` column covariance
#'   matrix.
#' @param lambda An optional numeric vector of length `p` containing finite,
#'   strictly positive exponential rate parameters. When `NULL`, a vector of
#'   ones is used.
#' @param log A logical scalar indicating whether the log-density should be
#'   returned. The default is `FALSE`.
#' @param eig_floor A positive numeric scalar specifying the minimum numerical
#'   eigenvalue allowed for the matrix `Q`. The default is `1e-10`.
#' @param q_policy Character string indicating how to handle cases in which
#'   the smallest eigenvalue of \code{Q} is less than \code{eig_floor}.
#'   Available options are \code{"warn"}, which regularizes \code{Q} and
#'   issues a warning; \code{"strict"}, which stops with an error; and
#'   \code{"regularize"}, which applies the regularization silently.
#'
#' @return A numeric scalar containing the MVREN density evaluated at `X`.
#'   When `log = TRUE`, the corresponding log-density is returned.
#'
#' @export
dmvren <- function(X,
                   M,
                   A,
                   Sigma,
                   Psi,
                   lambda = NULL,
                   log = FALSE,
                   eig_floor = 1e-10,
                   q_policy = c("warn", "strict", "regularize")) {
  q_policy <- match.arg(q_policy)
  X <- as.matrix(X)
  pars <- mvren_validate_parameters(
    X = X,
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda
  )

  p <- pars$p
  q <- pars$q
  Sigma_inv <- mvren_chol_solve(pars$Sigma, diag(p), "Sigma")
  Psi_inv <- mvren_chol_solve(pars$Psi, diag(q), "Psi")

  Y <- X - pars$M
  c_value <- sum(diag(Sigma_inv %*% Y %*% Psi_inv %*% t(Y)))
  b <- diag(pars$A %*% Psi_inv %*% t(Y) %*% Sigma_inv)
  Q <- Sigma_inv * (pars$A %*% Psi_inv %*% t(pars$A))
  Q <- mvren_prepare_Q(Q, policy = q_policy, eig_floor = eig_floor)$Q
  Q_inv <- mvren_chol_solve(Q, diag(p), "Q")

  h <- as.numeric(b - pars$lambda)
  mu <- as.numeric(Q_inv %*% h)

  log_density <-
    sum(log(pars$lambda)) -
    (p * (q - 1) / 2) * log(2 * pi) -
    (q / 2) * mvren_logdet(pars$Sigma, "Sigma") -
    (p / 2) * mvren_logdet(pars$Psi, "Psi") -
    (1 / 2) * mvren_logdet(Q, "Q") -
    (1 / 2) * c_value +
    (1 / 2) * as.numeric(crossprod(h, Q_inv %*% h)) +
    mvren_log_pmvnorm_upper(mu, mean = rep(0, p), sigma = Q_inv)

  if (isTRUE(log)) log_density else exp(log_density)
}
#' Observed-data log-likelihood of the MVREN distribution
#'
#' Evaluates the observed-data log-likelihood of a sample from the
#' matrix-variate row exponential-normal (MVREN) distribution by summing the
#' individual log-density contributions of the matrix-valued observations.
#'
#' @param X_array A numeric three-dimensional array with dimensions `p` by
#'   `q` by `n`, containing `n` matrix-valued observations.
#' @param M A numeric `p` by `q` location matrix.
#' @param A A numeric `p` by `q` row-specific latent-effect matrix.
#' @param Sigma A symmetric positive-definite `p` by `p` row covariance
#'   matrix.
#' @param Psi A symmetric positive-definite `q` by `q` column covariance
#'   matrix.
#' @param lambda An optional numeric vector of length `p` containing finite,
#'   strictly positive exponential rate parameters. When `NULL`, a vector of
#'   ones is used.
#' @param eig_floor A positive numeric scalar specifying the minimum numerical
#'   eigenvalue allowed for the matrix `Q` used in each density evaluation.
#'   The default is `1e-10`.
#' @param q_policy A character string specifying how to handle a matrix `Q`
#'   whose smallest eigenvalue is below `eig_floor`. Available options are
#'   `"warn"`, `"strict"`, and `"regularize"`; see [dmvren()] for details.
#'
#' @return A numeric scalar containing the observed-data log-likelihood of the
#'   supplied sample.
#'
#' @export
loglik_mvren <- function(X_array,
                         M,
                         A,
                         Sigma,
                         Psi,
                         lambda = NULL,
                         eig_floor = 1e-10,
                         q_policy = c("warn", "strict", "regularize")) {
  q_policy <- match.arg(q_policy)
  X_array <- mvren_validate_sample_array(X_array)
  n <- dim(X_array)[3]

  sum(vapply(seq_len(n), function(i) {
    dmvren(
      X = X_array[, , i],
      M = M,
      A = A,
      Sigma = Sigma,
      Psi = Psi,
      lambda = lambda,
      log = TRUE,
      eig_floor = eig_floor,
      q_policy = q_policy
    )
  }, numeric(1)))
}

# ---------------------------------------------------------------------------
# ECM estimation
# ---------------------------------------------------------------------------

#' Compute conditional expected residual cross-products
#'
#' Computes conditional expected row and column cross-product matrices for the
#' residual term `Y - WA`, using the conditional first and second moments of
#' the latent diagonal matrix `W`.
#'
#' @param Y A numeric `p` by `q` centered observation matrix.
#' @param A A numeric `p` by `q` latent-effect matrix.
#' @param wbar A numeric vector of length `p` containing the conditional means
#'   of the latent variables.
#' @param second_moment A numeric `p` by `p` matrix containing the conditional
#'   second moments of the latent variables.
#' @param Sigma_inv An optional numeric `p` by `p` inverse row covariance
#'   matrix. When supplied, the column cross-product matrix is computed.
#' @param Psi_inv An optional numeric `q` by `q` inverse column covariance
#'   matrix. When supplied, the row cross-product matrix is computed.
#'
#' @return A list containing any of the following components:
#' \describe{
#'   \item{column}{The conditional expectation of
#'   `t(Y - WA) %*% Sigma_inv %*% (Y - WA)`.}
#'   \item{row}{The conditional expectation of
#'   `(Y - WA) %*% Psi_inv %*% t(Y - WA)`.}
#' }
#' A component is included only when its corresponding inverse covariance
#' matrix is supplied.
#'
#' @keywords internal
#' @noRd
mvren_expected_crossproducts <- function(Y,
                                         A,
                                         wbar,
                                         second_moment,
                                         Sigma_inv = NULL,
                                         Psi_inv = NULL) {
  p <- nrow(Y)
  Wbar <- diag(wbar, p, p)
  output <- list()

  if (!is.null(Sigma_inv)) {
    # E[(Y - W A)' Sigma^{-1} (Y - W A) | X]
    output$column <-
      t(Y) %*% Sigma_inv %*% Y -
      t(Y) %*% Sigma_inv %*% Wbar %*% A -
      t(A) %*% Wbar %*% Sigma_inv %*% Y +
      t(A) %*% (Sigma_inv * second_moment) %*% A
    output$column <- (output$column + t(output$column)) / 2
  }

  if (!is.null(Psi_inv)) {
    # E[(Y - W A) Psi^{-1} (Y - W A)' | X]
    APsiA <- A %*% Psi_inv %*% t(A)
    output$row <-
      Y %*% Psi_inv %*% t(Y) -
      Y %*% Psi_inv %*% t(A) %*% Wbar -
      Wbar %*% A %*% Psi_inv %*% t(Y) +
      second_moment * APsiA
    output$row <- (output$row + t(output$row)) / 2
  }

  output
}
#' ECM estimation for the MVREN distribution
#'
#' Fits a matrix-variate row exponential-normal model by an expectation-
#' conditional-maximization algorithm with optional numerical regularization
#' and monotonicity-preserving step halving.
#'
#' @param X A numeric three-dimensional sample array with dimensions `p` by
#'   `q` by `n`.
#' @param max_iter A positive integer specifying the maximum number of ECM
#'   iterations. The default is `200`.
#' @param tol A positive numeric scalar specifying the relative log-likelihood
#'   convergence tolerance. The default is `1e-6`.
#' @param normalize_Psi A logical scalar indicating whether `Psi` should be
#'   rescaled toward a unit-determinant representation. Eigenvalue flooring
#'   applied after rescaling may produce a small numerical deviation from
#'   `det(Psi) = 1`. The default is `TRUE`.
#' @param lambda_mode A character string specifying the treatment of the rate
#'   vector. `"fixed"` keeps `lambda` fixed; `"unit_A_rows"` estimates
#'   `lambda` and normalizes every row of `A` to unit Euclidean norm; and
#'   `"unconstrained"` estimates `lambda` without a scale-identifying
#'   restriction and is therefore non-identifiable.
#' @param M_init An optional `p` by `q` starting value for `M`.
#' @param A_init An optional `p` by `q` starting value for `A`.
#' @param Sigma_init An optional positive-definite `p` by `p` starting value
#'   for `Sigma`.
#' @param Psi_init An optional positive-definite `q` by `q` starting value for
#'   `Psi`.
#' @param lambda_init An optional positive vector of length `p` containing the
#'   initial rate parameters.
#' @param verbose A logical scalar indicating whether iteration information
#'   should be printed. The default is `TRUE`.
#' @param eig_floor A positive numeric scalar used as an eigenvalue floor in
#'   covariance and precision-matrix regularization. The default is `1e-8`.
#' @param q_policy A character string specifying how to handle a matrix `Q`
#'   whose smallest eigenvalue is below `eig_floor`. Available options are
#'   `"warn"`, `"strict"`, and `"regularize"`.
#' @param monotone_tol A nonnegative numeric scalar specifying the largest
#'   numerical decrease in log-likelihood treated as acceptable. The default
#'   is `1e-7`.
#' @param max_backtracks A nonnegative integer specifying the maximum number of
#'   step halvings attempted when a full update decreases the log-likelihood.
#'   The default is `12L`.
#' @param moment_max_retries A positive integer specifying the maximum number
#'   of repeated attempts used to obtain valid truncated-normal moments. The
#'   default is `4L`.
#'
#' @return An object of class `"mvren_fit"`, represented by a list containing:
#' \describe{
#'   \item{M, A, Sigma, Psi, lambda}{The fitted MVREN parameters.}
#'   \item{loglik}{The final observed-data log-likelihood.}
#'   \item{loglik_history}{The sequence of accepted log-likelihood values.}
#'   \item{iterations}{The number of ECM iterations performed.}
#'   \item{BIC}{The Bayesian information criterion computed by the function.}
#'   \item{converged}{A logical indicator of convergence under `tol`.}
#'   \item{criterion}{The final relative log-likelihood change.}
#'   \item{monotone}{A logical indicator of whether decreases larger than
#'   `monotone_tol` were avoided.}
#'   \item{monotone_drops}{Any recorded decreases exceeding
#'   `monotone_tol`.}
#'   \item{step_sizes}{The accepted step size at each iteration.}
#'   \item{rejected_updates}{The number of full update sequences rejected
#'   after unsuccessful backtracking.}
#'   \item{normalize_Psi}{The requested covariance-identification setting.}
#'   \item{lambda_mode}{The selected rate-identification mode.}
#'   \item{npar}{The parameter count used in the BIC calculation.}
#' }
#'
#' @details
#' In the E-step, each latent vector is represented by a multivariate normal
#' distribution with covariance `Q^{-1}` and mean
#' `Q^{-1}(b_i - lambda)`, truncated to the positive orthant. The subsequent
#' CM-steps update `M`, `A`, `Psi`, `Sigma`, and, when requested, `lambda`.
#' Candidate updates that reduce the observed-data log-likelihood beyond
#' `monotone_tol` are subjected to step halving. Numerical eigenvalue flooring
#' can make the implemented procedure a stabilized ECM-type algorithm rather
#' than an exact unregularized ECM iteration.
#'
#' @keywords internal
mvren_ecm <- function(X,
                      max_iter = 200,
                      tol = 1e-6,
                      normalize_Psi = TRUE,
                      lambda_mode = c("fixed", "unit_A_rows", "unconstrained"),
                      M_init = NULL,
                      A_init = NULL,
                      Sigma_init = NULL,
                      Psi_init = NULL,
                      lambda_init = NULL,
                      verbose = TRUE,
                      eig_floor = 1e-8,
                      q_policy = c("warn", "strict", "regularize"),
                      monotone_tol = 1e-7,
                      max_backtracks = 12L,
                      moment_max_retries = 4L) {
  X <- mvren_validate_sample_array(X)
  lambda_mode <- match.arg(lambda_mode)
  q_policy <- match.arg(q_policy)
  max_backtracks <- as.integer(max_backtracks)
  if (length(max_backtracks) != 1L || !is.finite(max_backtracks) ||
      max_backtracks < 0L) {
    stop("max_backtracks must be a non-negative integer.", call. = FALSE)
  }

  if (lambda_mode == "unconstrained") {
    warning(
      paste0(
        "Estimating A and lambda without a scale constraint is non-identifiable. ",
        "Prefer lambda_mode = 'fixed' or 'unit_A_rows'."
      ),
      call. = FALSE
    )
  }

  p <- dim(X)[1]
  q <- dim(X)[2]
  n <- dim(X)[3]

  init <- mvren_initial_values(
    X = X,
    M_init = M_init,
    A_init = A_init,
    Sigma_init = Sigma_init,
    Psi_init = Psi_init,
    lambda_init = lambda_init,
    lambda_mode = lambda_mode
  )

  M <- init$M
  A <- init$A
  Sigma <- init$Sigma
  Psi <- init$Psi
  lambda <- init$lambda

  if (normalize_Psi) {
    scale_Psi <- det(Psi)^(1 / q)
    Psi <- Psi / scale_Psi
    Sigma <- Sigma * scale_Psi
  }

  current_loglik <- loglik_mvren(
    X,
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda,
    eig_floor = eig_floor,
    q_policy = q_policy
  )

  loglik_history <- current_loglik
  criterion <- Inf
  count <- 0L
  monotone <- TRUE
  monotone_drops <- numeric(0)
  step_sizes <- numeric(0)
  rejected_updates <- 0L
  update_rejected <- FALSE

  while (criterion > tol && count < max_iter) {
    count <- count + 1L

    Sigma <- mvren_make_posdef(Sigma, eig_floor = eig_floor)
    Psi <- mvren_make_posdef(Psi, eig_floor = eig_floor)
    M_old <- M
    A_old <- A
    Sigma_old <- Sigma
    Psi_old <- Psi
    lambda_old <- lambda
    Sigma_inv <- mvren_chol_solve(Sigma, diag(p), "Sigma")
    Psi_inv <- mvren_chol_solve(Psi, diag(q), "Psi")

    wbar_list <- vector("list", n)
    second_moment_list <- vector("list", n)

    # E-step: W_i | X_i follows a multivariate normal distribution with
    # mean Q^{-1}(b_i - lambda) and covariance Q^{-1}, truncated to R_+^p.
    Q <- Sigma_inv * (A %*% Psi_inv %*% t(A))
    Q <- mvren_prepare_Q(Q, policy = q_policy, eig_floor = eig_floor)$Q
    Q_inv <- mvren_chol_solve(Q, diag(p), "Q")

    for (i in seq_len(n)) {
      Y <- X[, , i] - M
      b <- diag(A %*% Psi_inv %*% t(Y) %*% Sigma_inv)
      conditional_mean <- as.numeric(Q_inv %*% (b - lambda))

      moments <- mvren_truncated_moments(
        mean = conditional_mean,
        sigma = Q_inv,
        eig_floor = eig_floor,
        max_retries = moment_max_retries
      )

      wbar_list[[i]] <- moments$mean
      second_moment_list[[i]] <- moments$second_moment
    }

    # CM-step 1: update M.
    M_new <- matrix(0, p, q)
    for (i in seq_len(n)) {
      Wbar <- diag(wbar_list[[i]], p, p)
      M_new <- M_new + X[, , i] - Wbar %*% A
    }
    M_new <- M_new / n

    # CM-step 2: update A.
    left <- matrix(0, p, p)
    right <- matrix(0, p, q)

    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      Wbar <- diag(wbar_list[[i]], p, p)
      left <- left + Sigma_inv * second_moment_list[[i]]
      right <- right + Wbar %*% Sigma_inv %*% Y
    }

    left <- mvren_make_posdef(left, eig_floor = eig_floor)
    A_new <- mvren_chol_solve(left, right, "A update")

    # CM-step 3: update Psi using the expected residual column crossproducts.
    # This expression is algebraically equivalent to the paper's Delta_i / B_ij
    # Cholesky representation, but avoids factorizing a potentially singular
    # positive-semidefinite pq x pq matrix.
    Psi_new <- matrix(0, q, q)
    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      crossproducts <- mvren_expected_crossproducts(
        Y = Y,
        A = A_new,
        wbar = wbar_list[[i]],
        second_moment = second_moment_list[[i]],
        Sigma_inv = Sigma_inv
      )
      Psi_new <- Psi_new + crossproducts$column
    }
    Psi_new <- mvren_make_posdef(Psi_new / (n * p), eig_floor = eig_floor)

    if (normalize_Psi) {
      scale_Psi <- det(Psi_new)^(1 / q)
      if (!is.finite(scale_Psi) || scale_Psi <= 0) {
        stop("The Psi update has a non-positive determinant.", call. = FALSE)
      }
      Psi_new <- mvren_make_posdef(Psi_new / scale_Psi, eig_floor = eig_floor)
    }

    # CM-step 4: update Sigma conditional on the updated Psi.
    Psi_new_inv <- mvren_chol_solve(Psi_new, diag(q), "Psi update")
    Sigma_new <- matrix(0, p, p)
    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      crossproducts <- mvren_expected_crossproducts(
        Y = Y,
        A = A_new,
        wbar = wbar_list[[i]],
        second_moment = second_moment_list[[i]],
        Psi_inv = Psi_new_inv
      )
      Sigma_new <- Sigma_new + crossproducts$row
    }
    Sigma_new <- mvren_make_posdef(Sigma_new / (n * q), eig_floor = eig_floor)

    # CM-step 5: update lambda when requested.
    if (lambda_mode == "fixed") {
      lambda_new <- lambda
    } else {
      sum_wbar <- Reduce(`+`, wbar_list)
      if (any(!is.finite(sum_wbar)) || any(sum_wbar <= 0)) {
        stop("The lambda CM-step received invalid conditional latent means.", call. = FALSE)
      }
      lambda_new <- n / sum_wbar
    }

    # Apply the row-norm identifying restriction only after all CM updates.
    identified <- mvren_identify_A_lambda(
      A = A_new,
      lambda = lambda_new,
      mode = lambda_mode
    )
    A_new <- identified$A
    lambda_new <- identified$lambda

    evaluate_candidate <- function(M_value, A_value, Sigma_value, Psi_value,
                                   lambda_value) {
      moments_are_available <- tryCatch({
        Sigma_value_inv <- mvren_chol_solve(
          Sigma_value, diag(p), "candidate Sigma"
        )
        Psi_value_inv <- mvren_chol_solve(
          Psi_value, diag(q), "candidate Psi"
        )
        Q_value <- Sigma_value_inv *
          (A_value %*% Psi_value_inv %*% t(A_value))
        Q_value <- mvren_prepare_Q(
          Q_value, policy = q_policy, eig_floor = eig_floor
        )$Q
        Q_value_inv <- mvren_chol_solve(
          Q_value, diag(p), "candidate Q"
        )
        for (index in seq_len(n)) {
          Y_value <- X[, , index] - M_value
          b_value <- diag(
            A_value %*% Psi_value_inv %*% t(Y_value) %*% Sigma_value_inv
          )
          candidate_moments <- mvren_truncated_moments(
            mean = as.numeric(Q_value_inv %*% (b_value - lambda_value)),
            sigma = Q_value_inv, eig_floor = eig_floor,
            max_retries = moment_max_retries
          )
          if (any(!is.finite(candidate_moments$mean)) ||
              any(!is.finite(candidate_moments$variance))) stop("invalid moments")
        }
        TRUE
      }, error = function(e) FALSE)
      if (!moments_are_available) return(NA_real_)
      tryCatch(loglik_mvren(
        X, M = M_value, A = A_value, Sigma = Sigma_value,
        Psi = Psi_value, lambda = lambda_value, eig_floor = eig_floor,
        q_policy = q_policy
      ), error = function(e) NA_real_)
    }
    new_loglik <- evaluate_candidate(M_new, A_new, Sigma_new, Psi_new, lambda_new)
    step_size <- 1

    if (!is.finite(new_loglik) || new_loglik < current_loglik - monotone_tol) {
      accepted <- FALSE
      for (backtrack in seq_len(as.integer(max_backtracks))) {
        step_size <- 2^(-backtrack)
        M_try <- M_old + step_size * (M_new - M_old)
        A_try <- A_old + step_size * (A_new - A_old)
        Sigma_try <- mvren_make_posdef(
          Sigma_old + step_size * (Sigma_new - Sigma_old), eig_floor
        )
        Psi_try <- mvren_make_posdef(
          Psi_old + step_size * (Psi_new - Psi_old), eig_floor
        )
        if (normalize_Psi) Psi_try <- Psi_try / det(Psi_try)^(1 / q)
        lambda_try <- lambda_old + step_size * (lambda_new - lambda_old)
        identified_try <- mvren_identify_A_lambda(
          A_try, lambda_try, mode = lambda_mode
        )
        A_try <- identified_try$A
        lambda_try <- identified_try$lambda
        loglik_try <- evaluate_candidate(
          M_try, A_try, Sigma_try, Psi_try, lambda_try
        )
        if (is.finite(loglik_try) && loglik_try >= current_loglik - monotone_tol) {
          M_new <- M_try
          A_new <- A_try
          Sigma_new <- Sigma_try
          Psi_new <- Psi_try
          lambda_new <- lambda_try
          new_loglik <- loglik_try
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        M_new <- M_old
        A_new <- A_old
        Sigma_new <- Sigma_old
        Psi_new <- Psi_old
        lambda_new <- lambda_old
        new_loglik <- current_loglik
        step_size <- 0
        rejected_updates <- rejected_updates + 1L
        update_rejected <- TRUE
      }
    }

    difference <- new_loglik - current_loglik
    if (difference < -monotone_tol) {
      monotone <- FALSE
      monotone_drops <- c(monotone_drops, difference)
    }
    step_sizes <- c(step_sizes, step_size)

    criterion <- abs(difference) / (abs(current_loglik) + .Machine$double.eps)

    M <- M_new
    A <- A_new
    Sigma <- Sigma_new
    Psi <- Psi_new
    lambda <- lambda_new
    current_loglik <- new_loglik
    loglik_history <- c(loglik_history, current_loglik)

    if (verbose) {
      cat(sprintf(
        "Iteration: %d | logLik: %.8f | relative change: %.3e\n",
        count,
        current_loglik,
        criterion
      ))
    }
    if (update_rejected) break
  }

  converged <- criterion <= tol && !update_rejected
  if (!converged && count == max_iter) {
    warning(
      "The ECM algorithm reached max_iter before satisfying the convergence criterion.",
      call. = FALSE
    )
  }

  covariance_parameters <-
    p * (p + 1) / 2 + q * (q + 1) / 2 - as.integer(normalize_Psi)

  if (lambda_mode == "fixed") {
    npar <- 2 * p * q + covariance_parameters
  } else if (lambda_mode == "unit_A_rows") {
    # p row-norm restrictions on A offset the p estimated rate parameters.
    npar <- 2 * p * q + covariance_parameters
  } else {
    npar <- 2 * p * q + covariance_parameters + p
  }

  BIC <- -2 * current_loglik + npar * log(n)


  obj.out <- list(M = M,
                  A = A,
                  Sigma = Sigma,
                  Psi = Psi,
                  lambda = lambda,
                  loglik = current_loglik,
                  loglik_history = loglik_history,
                  iterations = count,
                  BIC = BIC,
                  converged = converged,
                  criterion = criterion,
                  monotone = monotone,
                  monotone_drops = monotone_drops,
                  step_sizes = step_sizes,
                  rejected_updates = rejected_updates,
                  normalize_Psi = normalize_Psi,
                  lambda_mode = lambda_mode,
                  npar = npar)

    class(obj.out) <- "MVREN.ECM"
    obj.out
}

#' Print an MVREN model fit
#'
#' Prints a concise summary of an MVREN ECM fit. This function is retained as
#' an internal helper and is not documented or exported as part of the public
#' package interface.
#'
#' @param x An object of class `"mvren_fit"`.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @keywords internal
#' @noRd
print.mvren_fit <- function(x, ...) {
  cat("Matrix-Variate Row Exponential-Normal ECM fit\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = 8), "\n")
  cat("BIC:", format(x$BIC, digits = 8), "\n")
  cat("Lambda mode:", x$lambda_mode, "\n")
  cat("lambda:", paste(format(x$lambda, digits = 6), collapse = ", "), "\n")
  invisible(x)
}

#' Predefined parameter set for the MVREN distribution
#'
#' Constructs a predefined collection of parameters for the matrix-variate
#' row exponential-normal distribution, intended for simulation studies,
#' numerical experiments, and reproducible examples.
#'
#' @param normalize_Psi Logical scalar indicating whether \code{Psi} should be
#'   normalized to have determinant one. The corresponding reciprocal scaling
#'   is applied to \code{Sigma} so that the Kronecker-product covariance
#'   structure remains unchanged. The default is \code{TRUE}.
#'
#' @return A list containing the matrices \code{M}, \code{A}, \code{Sigma},
#'   and \code{Psi}, together with the exponential-rate vector \code{lambda}.
#'
#' @export
mvren_article_parameters <- function(normalize_Psi = TRUE) {
  M <- matrix(c(
    0.50, 1.00, 1.00, 1.00,
    1.50, 1.00, 1.50, 0.50,
    0.50, 1.00, 0.50, 1.50
  ), nrow = 3, ncol = 4, byrow = TRUE)

  A <- matrix(c(
    0.18, 0.25, 0.71, 1.62,
    0.73, 1.42, 1.16, 0.83,
    0.22, 0.02, 0.13, 0.81
  ), nrow = 3, ncol = 4, byrow = TRUE)

  Sigma <- matrix(c(
    0.100, 0.040, 0.016,
    0.040, 0.100, 0.040,
    0.016, 0.040, 0.100
  ), nrow = 3, ncol = 3, byrow = TRUE)

  Psi <- matrix(c(
    2.152, 1.721, 1.377, 1.102,
    1.721, 2.152, 1.721, 1.377,
    1.377, 1.721, 2.152, 1.721,
    1.102, 1.377, 1.721, 2.152
  ), nrow = 4, ncol = 4, byrow = TRUE)

  if (normalize_Psi) {
    scale_Psi <- det(Psi)^(1 / nrow(Psi))
    Psi <- Psi / scale_Psi
    Sigma <- Sigma * scale_Psi
  }

  list(
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = rep(1, nrow(M))
  )
}

#' Compute MVREN fit diagnostics
#'
#' Computes estimation-error, convergence, likelihood, and identifiability
#' diagnostics for a fitted MVREN model relative to a supplied set of true
#' parameter values.
#'
#' @param fit A fitted MVREN model object containing the estimated parameters
#'   `M`, `A`, `Sigma`, `Psi`, and `lambda`, together with the components
#'   `loglik`, `BIC`, `iterations`, `converged`, and `monotone`.
#' @param truth A list containing the true parameter values `M`, `A`, `Sigma`,
#'   `Psi`, and `lambda`.
#'
#' @return A one-row data frame containing:
#' \describe{
#'   \item{err_M}{Relative Frobenius error for the location matrix.}
#'   \item{err_A}{Relative Frobenius error for the latent-effect matrix.}
#'   \item{err_Sigma}{Relative Frobenius error for the row covariance matrix.}
#'   \item{err_Psi}{Relative Frobenius error for the column covariance matrix.}
#'   \item{err_lambda}{Relative Euclidean error for the rate vector.}
#'   \item{loglik}{Observed-data log-likelihood of the fitted model.}
#'   \item{BIC}{Bayesian information criterion of the fitted model.}
#'   \item{iterations}{Number of iterations performed by the fitting
#'   algorithm.}
#'   \item{converged}{Logical indicator of algorithmic convergence.}
#'   \item{monotone}{Logical indicator of monotonic log-likelihood behavior.}
#'   \item{det_Psi}{Determinant of the estimated column covariance matrix.}
#' }
#'
#' @keywords internal
#' @noRd
mvren_fit_diagnostics <- function(fit, truth) {
  data.frame(
    err_M = mvren_relative_frobenius(fit$M, truth$M),
    err_A = mvren_relative_frobenius(fit$A, truth$A),
    err_Sigma = mvren_relative_frobenius(fit$Sigma, truth$Sigma),
    err_Psi = mvren_relative_frobenius(fit$Psi, truth$Psi),
    err_lambda = sqrt(sum((fit$lambda - truth$lambda)^2)) /
      max(sqrt(sum(truth$lambda^2)), .Machine$double.eps),
    loglik = fit$loglik,
    BIC = fit$BIC,
    iterations = fit$iterations,
    converged = fit$converged,
    monotone = fit$monotone,
    det_Psi = det(fit$Psi)
  )
}
#' Run a reproducible Monte Carlo audit of the MVREN ECM estimator
#'
#' Performs a Monte Carlo simulation study to evaluate the finite-sample
#' behavior of the MVREN ECM estimator across a collection of sample sizes.
#' For each replication, data are generated from a specified MVREN parameter
#' configuration, the model is fitted, and estimation, convergence, and
#' likelihood diagnostics are recorded.
#'
#' @param sample_sizes A numeric vector of positive whole-number sample sizes
#'   to be evaluated. The default is `c(50, 100, 200, 400, 800, 1600)`.
#' @param replications A positive whole number specifying the number of Monte
#'   Carlo replications performed for each sample size. The default is `200`.
#' @param truth A list containing the data-generating MVREN parameters `M`,
#'   `A`, `Sigma`, `Psi`, and `lambda`. By default, the parameters returned by
#'   `mvren_article_parameters()` are used.
#' @param lambda_mode A character string specifying the identifiability
#'   convention used when fitting the model. Supported values are `"fixed"`,
#'   `"unit_A_rows"`, and `"unconstrained"`. The default is `"fixed"`.
#' @param max_iter A positive integer specifying the maximum number of ECM
#'   iterations allowed for each fitted model. The default is `200`.
#' @param tol A positive numeric scalar specifying the convergence tolerance
#'   passed to the ECM algorithm. The default is `1e-6`.
#' @param seed An integer used to initialize the random-number generator,
#'   ensuring reproducibility. The default is `123`.
#' @param verbose A logical scalar indicating whether simulation progress and
#'   ECM output should be printed. The default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{A data frame with one row per replication, containing the
#'   sample size, replication number, relative estimation errors, likelihood
#'   diagnostics, convergence information, determinant of the estimated
#'   `Psi`, and any fitting error message.}
#'   \item{summary}{A data frame summarizing the results for each sample size,
#'   including convergence and monotonicity rates, median relative estimation
#'   errors, and median iteration counts.}
#'   \item{generating_truth}{The parameter set used to generate the simulated
#'   observations.}
#'   \item{comparison_truth}{The parameter set used to compute estimation
#'   errors. Under `lambda_mode = "unit_A_rows"`, this contains the
#'   identifiability-adjusted values of `A` and `lambda`.}
#' }
#'
#' @export
mvren_monte_carlo <- function(sample_sizes = c(50, 100, 200, 400, 800, 1600),
                              replications = 200,
                              truth = mvren_article_parameters(),
                              lambda_mode = "fixed",
                              max_iter = 200,
                              tol = 1e-6,
                              seed = 123,
                              verbose = FALSE) {
  if (any(sample_sizes <= 0L) || replications <= 0L) {
    stop("sample_sizes and replications must be positive.", call. = FALSE)
  }

  mvren_validate_parameters(
    M = truth$M,
    A = truth$A,
    Sigma = truth$Sigma,
    Psi = truth$Psi,
    lambda = truth$lambda
  )

  comparison_truth <- truth
  if (lambda_mode == "unit_A_rows") {
    identified_truth <- mvren_identify_A_lambda(
      A = truth$A,
      lambda = truth$lambda,
      mode = "unit_A_rows"
    )
    comparison_truth$A <- identified_truth$A
    comparison_truth$lambda <- identified_truth$lambda
  }

  set.seed(seed)
  rows <- vector("list", length(sample_sizes) * replications)
  index <- 0L

  for (n in sample_sizes) {
    for (replication in seq_len(replications)) {
      index <- index + 1L

      if (verbose) {
        cat(sprintf(
          "[MVREN Monte Carlo] n = %d, replication = %d/%d\n",
          n,
          replication,
          replications
        ))
      }

      X <- rmvren(
        n = n,
        M = truth$M,
        A = truth$A,
        Sigma = truth$Sigma,
        Psi = truth$Psi,
        lambda = truth$lambda
      )

      fit <- tryCatch(
        mvren_ecm(
          X = X,
          max_iter = max_iter,
          tol = tol,
          lambda_mode = lambda_mode,
          lambda_init = truth$lambda,
          verbose = verbose
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        rows[[index]] <- data.frame(
          n = n,
          replication = replication,
          err_M = NA_real_,
          err_A = NA_real_,
          err_Sigma = NA_real_,
          err_Psi = NA_real_,
          err_lambda = NA_real_,
          loglik = NA_real_,
          BIC = NA_real_,
          iterations = NA_integer_,
          converged = FALSE,
          monotone = FALSE,
          det_Psi = NA_real_,
          error = conditionMessage(fit)
        )
      } else {
        diagnostics <- mvren_fit_diagnostics(fit, comparison_truth)
        rows[[index]] <- cbind(
          data.frame(n = n, replication = replication),
          diagnostics,
          data.frame(error = NA_character_)
        )
      }
    }
  }

  results <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(results, results$n), function(df) {
    data.frame(
      n = unique(df$n),
      replications = nrow(df),
      convergence_rate = mean(df$converged, na.rm = TRUE),
      monotone_rate = mean(df$monotone, na.rm = TRUE),
      median_err_M = stats::median(df$err_M, na.rm = TRUE),
      median_err_A = stats::median(df$err_A, na.rm = TRUE),
      median_err_Sigma = stats::median(df$err_Sigma, na.rm = TRUE),
      median_err_Psi = stats::median(df$err_Psi, na.rm = TRUE),
      median_err_lambda = stats::median(df$err_lambda, na.rm = TRUE),
      median_iterations = stats::median(df$iterations, na.rm = TRUE)
    )
  }))
  rownames(summary) <- NULL

  list(
    results = results,
    summary = summary,
    generating_truth = truth,
    comparison_truth = comparison_truth
  )
}


# ---------------------------------------------------------------------------
# Standalone example execution block
#
# The top-level calls below run when this file is sourced. The call to
# print(fit) requires an existing object named fit; remove or comment out
# these calls before using this file as package source code.
# ---------------------------------------------------------------------------

mvren_run_example <- function(seed = 123L, n = 100L, max_iter = 200L,
                              verbose = TRUE) {
  set.seed(seed)
  pars <- mvren_article_parameters()
  X <- rmvren(n = n, M = pars$M, A = pars$A, Sigma = pars$Sigma,
              Psi = pars$Psi, lambda = pars$lambda)
  fit <- mvren_ecm(
    X, max_iter = max_iter, lambda_mode = "fixed",
    lambda_init = pars$lambda, verbose = verbose
  )
  print(fit)
  invisible(fit)
}
