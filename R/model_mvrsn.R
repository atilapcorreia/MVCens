# Matrix-Variate Row Skew-Normal (MVRSN)
#
# This file is intentionally self-contained. It implements the theoretical
# construction used in the paper:
#
#   X = M + W A + V,
#
# where W = diag(W_1, ..., W_p), W_i are independent standard half-normal
# variables, and V is matrix-normal with row covariance Sigma and column
# covariance Psi. The public functions below cover the full empirical loop:
# density -> likelihood -> simulation -> ECM estimation -> Monte Carlo audit.


# ---------------------------------------------------------------------------
# Internal numerical utilities
# ---------------------------------------------------------------------------

#' Compute the lower Cholesky factor
#'
#' Internal utility that computes the lower triangular Cholesky factor of a
#' symmetric positive definite matrix. An error is raised if the input matrix is
#' not positive definite.
#'
#' @param S A numeric symmetric positive definite matrix.
#' @param name A character string identifying the matrix in error messages.
#'   Defaults to `"matrix"`.
#'
#' @return A lower triangular matrix `L` such that `S = L %*% t(L)`.
#'
#' @keywords internal
#' @noRd
mvrsn_chol_lower <- function(S, name = "matrix") {
  S <- as.matrix(S)
  matrix_assert_posdef(S, name)
  t(chol(S))
}

#' Validate model parameters
#'
#' Internal utility that validates the dimensions and numeric types of the
#' model parameters. Optionally, the covariance matrices are checked for
#' symmetric positive definiteness. The validated objects are returned in a
#' standardized format for subsequent computations.
#'
#' @param X An optional numeric matrix to be checked against the dimensions of
#'   `M` and `A`.
#' @param M A numeric location matrix.
#' @param A A numeric skewness matrix.
#' @param Sigma A numeric row covariance matrix.
#' @param Psi A numeric column covariance matrix.
#' @param check_posdef Logical indicating whether `Sigma` and `Psi` should be
#'   checked for symmetric positive definiteness. Defaults to `TRUE`.
#'
#' @return A list containing the validated matrices `M`, `A`, `Sigma`, and
#'   `Psi`, together with the corresponding matrix dimensions `p` and `q`.
#'
#' @keywords internal
#' @noRd
mvrsn_validate_parameters <- function(X = NULL, M, A, Sigma, Psi, check_posdef = TRUE) {
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
    stop("A must have the same dimension as M.", call. = FALSE)
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
      stop("X must have the same dimension as M and A.", call. = FALSE)
    }
  }

  if (check_posdef) {
    matrix_assert_posdef(Sigma, "Sigma")
    matrix_assert_posdef(Psi, "Psi")
  }

  list(M = M, A = A, Sigma = Sigma, Psi = Psi, p = p, q = q)
}

#' Validate a sample array
#'
#' Internal utility that verifies whether the input is a valid three-dimensional
#' array representing a sample of matrix-valued observations. The array is
#' required to have dimensions `p × q × n`, where `p` and `q` denote the matrix
#' dimensions and `n` is the sample size.
#'
#' @param X An object representing a sample of matrix-valued observations.
#'
#' @return A three-dimensional array with dimensions `p × q × n`.
#'
#' @keywords internal
#' @noRd
mvrsn_validate_sample_array <- function(X) {
  X <- as.array(X)
  dims <- dim(X)
  if (length(dims) != 3L) {
    stop("X must be a 3-dimensional array with dimensions p x q x n.", call. = FALSE)
  }
  if (any(dims <= 0L)) {
    stop("X dimensions must all be positive.", call. = FALSE)
  }
  X
}

#' Construct the diagonal selection matrix
#'
#' Internal utility that constructs the matrix used to extract or represent the
#' diagonal entries of a square matrix in vectorized form. The returned matrix
#' has as its columns the vectorized canonical diagonal basis matrices.
#'
#' @param p A positive integer specifying the dimension of the square matrices.
#'
#' @return A matrix of dimension `p^2 × p` whose `i`-th column is the
#'   vectorization of the `p × p` matrix having a one in the `(i,i)` position
#'   and zeros elsewhere.
#'
#' @keywords internal
#' @noRd
mvrsn_make_F <- function(p) {
  F <- matrix(0, p * p, p)
  for (i in seq_len(p)) {
    E <- matrix(0, p, p)
    E[i, i] <- 1
    F[, i] <- matrix_vectorize(E)
  }
  F
}

#' Generate initial parameter values
#'
#' Internal utility that constructs initial values for the ECM algorithm.
#' Unless supplied by the user, the location matrix is initialized from the
#' sample mean, the skewness matrix is initialized from standardized sample
#' third moments, and the row and column covariance matrices are initialized as
#' identity matrices. A small deterministic perturbation is used whenever the
#' estimated skewness matrix is numerically zero to avoid a degenerate starting
#' point.
#'
#' @param X A three-dimensional array of matrix-valued observations with
#'   dimensions `p × q × n`.
#' @param M_init An optional initial value for the location matrix.
#' @param A_init An optional initial value for the skewness matrix.
#' @param Sigma_init An optional initial value for the row covariance matrix.
#' @param Psi_init An optional initial value for the column covariance matrix.
#' @param skew_scale A positive numeric constant controlling the scale of the
#'   automatically generated initial skewness matrix. Defaults to `0.1`.
#'
#' @return A list containing the initial parameter values `M`, `A`, `Sigma`,
#'   and `Psi`.
#'
#' @keywords internal
#' @noRd
mvrsn_initial_values <- function(X,
                                 M_init = NULL,
                                 A_init = NULL,
                                 Sigma_init = NULL,
                                 Psi_init = NULL,
                                 skew_scale = 0.1) {
  X <- mvrsn_validate_sample_array(X)
  p <- dim(X)[1]
  q <- dim(X)[2]

  Xbar <- apply(X, c(1, 2), mean)
  M <- if (is.null(M_init)) Xbar else as.matrix(M_init)

  if (is.null(A_init)) {
    centered <- X
    for (i in seq_len(dim(X)[3])) {
      centered[, , i] <- X[, , i] - Xbar
    }

    sds <- apply(X, c(1, 2), stats::sd)
    sds[!is.finite(sds) | sds <= .Machine$double.eps] <- 1
    skew <- apply(centered, c(1, 2), function(z) mean(z^3)) / (sds^3)
    skew[!is.finite(skew)] <- 0

    A <- skew_scale * skew

    # Symmetric or standardized data can have sample skewness very close to
    # zero. A numerically zero A is an absorbing point of the ECM update, so we
    # use a tiny deterministic pattern as a last-resort perturbation.
    if (frobenius_norm(A) <= .Machine$double.eps^0.25) {
      pattern <- matrix(seq_len(p * q), nrow = p, ncol = q)
      pattern <- pattern - mean(pattern)
      A <- skew_scale * pattern / frobenius_norm(pattern)
    }

    M <- Xbar - sqrt(2 / pi) * A
  } else {
    A <- as.matrix(A_init)
  }

  Sigma <- if (is.null(Sigma_init)) diag(p) else as.matrix(Sigma_init)
  Psi <- if (is.null(Psi_init)) diag(q) else as.matrix(Psi_init)

  list(M = M, A = A, Sigma = Sigma, Psi = Psi)
}

#' Compute the log-probability of a multivariate normal distribution
#'
#' Internal utility that computes the logarithm of the cumulative probability of
#' a multivariate normal distribution over the region
#' \eqn{(-\infty, \mathrm{upper}]}. If the computed probability is numerically
#' zero or non-finite, the logarithm of the smallest positive normalized
#' floating-point number is returned to avoid numerical underflow.
#'
#' @param upper A numeric vector specifying the upper integration limits.
#' @param mean A numeric vector specifying the mean of the multivariate normal
#'   distribution.
#' @param sigma A symmetric positive definite covariance matrix.
#'
#' @return The logarithm of the corresponding multivariate normal cumulative
#'   probability.
#'
#' @keywords internal
#' @noRd
mvrsn_log_pmvnorm_upper <- function(upper, mean, sigma) {
  prob <- as.numeric(mvtnorm::pmvnorm(upper = upper, mean = mean, sigma = sigma))
  if (!is.finite(prob) || prob <= 0) {
    return(log(.Machine$double.xmin))
  }
  log(prob)
}


# ---------------------------------------------------------------------------
# Theoretical moments implied by X = M + W A + V
# ---------------------------------------------------------------------------

#' Mean matrix of the MVRSN distribution
#'
#' Computes the theoretical mean under the matrix-variate row skew-normal
#' (MVRSN) representation. In this model, independent standard half-normal
#' latent variables act row by row through a diagonal skewing mechanism,
#' producing the stochastic form `X = M + W A + V`, with matrix-normal noise
#' `V`.
#'
#' @param M Numeric `p` by `q` location matrix.
#' @param A Numeric `p` by `q` skewness matrix.
#'
#' @return A numeric `p` by `q` theoretical mean matrix.
#'
#' @examples
#' M <- matrix(0, 3, 4)
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' mvrsn_mean(M, A)
#' @family MVCens MVRSN functions
#' @export
mvrsn_mean <- function(M, A) {
  M <- as.matrix(M)
  A <- as.matrix(A)

  if (!all(dim(M) == dim(A))) {
    stop("M and A must have the same dimension.", call. = FALSE)
  }
  M + sqrt(2 / pi) * A
}

#' Row and column covariance summaries of the MVRSN distribution
#'
#' Computes the theoretical row and column covariance summaries of the
#' matrix-variate row skew-normal distribution. The summaries combine the
#' separable matrix-normal covariance contribution with the row-specific
#' half-normal skewing contribution induced by `A`. They correspond to row-wise
#' and column-wise second moments of the matrix-valued observation, not to the
#' full vectorized Kronecker covariance matrix.
#'
#' @param M Location matrix. It is accepted for API symmetry and dimension
#'   validation, but it does not enter the covariance formula.
#' @param A Skewness matrix of dimension \eqn{p \times q}.
#' @param Sigma Positive definite row covariance matrix of dimension
#'   \eqn{p \times p}.
#' @param Psi Positive definite column covariance matrix of dimension
#'   \eqn{q \times q}.
#' @return A list with `row` (`p` by `p`) and `column` (`q` by `q`) covariance
#'   summaries.
#'
#' @examples
#' M <- matrix(0, 3, 4)
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' mvrsn_covariances(
#'   M = M, A = A,
#'   Sigma = diag(3), Psi = diag(4)
#' )
#' @family MVCens MVRSN functions
#' @export
mvrsn_covariances <- function(M, A, Sigma, Psi) {
  pars <- mvrsn_validate_parameters(M = M, A = A, Sigma = Sigma, Psi = Psi)
  AAt <- pars$A %*% t(pars$A)
  list(
    row = (1 - 2 / pi) * diag(diag(AAt), pars$p, pars$p) + sum(diag(pars$Psi)) * pars$Sigma,
    column = (1 - 2 / pi) * t(pars$A) %*% pars$A + sum(diag(pars$Sigma)) * pars$Psi
  )
}

# ---------------------------------------------------------------------------
# Density and likelihood
# ---------------------------------------------------------------------------

#' Density of the matrix-variate row skew-normal distribution
#'
#' Evaluates the closed-form probability density or log-density of a
#' matrix-variate row skew-normal distribution at one `p` by `q`
#' matrix-valued observation. Unlike the common-latent MVSN model, MVRSN assigns
#' an independent half-normal skewing variable to each row, allowing
#' heterogeneous row-level asymmetry while retaining separable matrix-normal row
#' and column covariance components.
#'
#' @param X Matrix of dimension \eqn{p \times q}.
#' @param M Location matrix of dimension \eqn{p \times q}.
#' @param A Skewness matrix of dimension \eqn{p \times q}.
#' @param Sigma Positive definite row covariance matrix of dimension
#'   \eqn{p \times p}.
#' @param Psi Positive definite column covariance matrix of dimension
#'   \eqn{q \times q}.
#' @param log Logical. If \code{TRUE}, returns the log-density.
#'
#' @return A numeric scalar containing the density evaluated at \code{X}.
#'   When \code{log = TRUE}, the corresponding log-density is returned.
#'
#' @examples
#' X <- matrix(seq(-0.2, 0.9, length.out = 12), 3, 4)
#' M <- matrix(0, 3, 4)
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' dmvrsn(
#'   X = X, M = M, A = A,
#'   Sigma = diag(3), Psi = diag(4)
#' )
#' @family MVCens density functions
#' @family MVCens MVRSN functions
#' @export
dmvrsn <- function(X, M, A, Sigma, Psi, log = FALSE) {

  X <- as.matrix(X)
  pars <- mvrsn_validate_parameters(X = X, M = M, A = A, Sigma = Sigma, Psi = Psi)
  p <- pars$p
  q <- pars$q

  Sigma_inv <- solve(pars$Sigma)
  Psi_inv <- solve(pars$Psi)

  Y <- X - pars$M
  c_val <- sum(diag(Sigma_inv %*% Y %*% Psi_inv %*% t(Y)))
  b <- diag(pars$A %*% Psi_inv %*% t(Y) %*% Sigma_inv)
  Q <- diag(p) + Sigma_inv * (pars$A %*% Psi_inv %*% t(pars$A))

  matrix_assert_posdef(Q, "Q")
  Q_inv <- solve(Q)
  mu_cdf <- as.vector(Q_inv %*% b)

  log_det_Sigma <- as.numeric(determinant(pars$Sigma, logarithm = TRUE)$modulus)
  log_det_Psi <- as.numeric(determinant(pars$Psi, logarithm = TRUE)$modulus)
  log_det_Q <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  log_cdf <- mvrsn_log_pmvnorm_upper(mu_cdf, mean = rep(0, p), sigma = Q_inv)

  log_density <-
    p * log(2) -
    (p * q / 2) * log(2 * pi) -
    (q / 2) * log_det_Sigma -
    (p / 2) * log_det_Psi -
    (1 / 2) * log_det_Q -
    (1 / 2) * c_val +
    (1 / 2) * as.numeric(t(b) %*% Q_inv %*% b) +
    log_cdf

  if (log) log_density else exp(log_density)
}

#' Log-likelihood function of the MVRSN distribution
#'
#' Sums the closed-form MVRSN log-density over a sample array with dimensions
#' `p` by `q` by `n`. Each slice `X_array[, , i]` is treated as one
#' matrix-valued observation with row-specific latent skewness.
#'
#' @param X_array Numeric sample array with dimensions `p` by `q` by `n`.
#' @param M Numeric `p` by `q` location matrix.
#' @param A Numeric `p` by `q` skewness matrix.
#' @param Sigma Positive-definite `p` by `p` row covariance matrix.
#' @param Psi Positive-definite `q` by `q` column covariance matrix.
#' @return A finite numeric scalar containing the observed-data log-likelihood.
#'
#' @examples
#' x <- array(seq(-0.2, 2.1, length.out = 24), dim = c(3, 4, 2))
#' M <- matrix(0, 3, 4)
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' loglik_mvrsn(
#'   x, M = M, A = A,
#'   Sigma = diag(3), Psi = diag(4)
#' )
#' @family MVCens likelihood functions
#' @family MVCens MVRSN functions
#' @export
loglik_mvrsn <- function(X_array, M, A, Sigma, Psi) {
  X_array <- mvrsn_validate_sample_array(X_array)
  n <- dim(X_array)[3]

  sum(vapply(seq_len(n), function(i) {
    dmvrsn(X_array[, , i], M, A, Sigma, Psi, log = TRUE)
  }, numeric(1)))
}


# ---------------------------------------------------------------------------
# ECM estimation
# ---------------------------------------------------------------------------

#' ECM estimation for the MVRSN distribution
#'
#' Fits the MVRSN model by Expectation/Conditional Maximization. The E-step
#' uses the conditional truncated normal law of the latent vector W | X. The
#' CM-steps implement the updates derived in the paper for M, A, Psi, and
#' Sigma. The function records the likelihood path so the simulation layer can
#' check whether the numerical implementation behaves like an EM/ECM algorithm.
#'
#' @param X Sample array p x q x n.
#' @param max_iter Maximum ECM iterations.
#' @param precision Relative log-likelihood convergence tolerance.
#' @param normalize_Psi Enforce the identifiability constraint det(Psi) = 1.
#' @param M_init Optional initial M.
#' @param A_init Optional initial A.
#' @param Sigma_init Optional initial Sigma.
#' @param Psi_init Optional initial Psi.
#' @param verbose Print iteration log when TRUE.
#' @param eig_floor Minimum eigenvalue used in numerical positive-definite
#'   projections.
#' @param monotone_tol Allowed numerical drop in log-likelihood before a
#'   monotonicity flag is raised.
#' @return A list with estimates, log-likelihood history, convergence flags,
#'   BIC, and monotonicity diagnostics.
#'
#' @keywords internal
mvrsn_ecm <- function(X,
                      max_iter = 200,
                      precision = 1e-6,
                      normalize_Psi = TRUE,
                      M_init = NULL,
                      A_init = NULL,
                      Sigma_init = NULL,
                      Psi_init = NULL,
                      verbose = TRUE,
                      eig_floor = 1e-8,
                      monotone_tol = 1e-7) {

  X <- mvrsn_validate_sample_array(X)
  p <- dim(X)[1]
  q <- dim(X)[2]
  n <- dim(X)[3]

  init <- mvrsn_initial_values(
    X,
    M_init = M_init,
    A_init = A_init,
    Sigma_init = Sigma_init,
    Psi_init = Psi_init
  )
  M <- init$M
  A <- init$A
  Sigma <- init$Sigma
  Psi <- init$Psi

  mvrsn_validate_parameters(M = M, A = A, Sigma = Sigma, Psi = Psi)
  if (normalize_Psi) {
    Psi <- Psi / det(Psi)^(1 / q)
  }

  Fmat <- mvrsn_make_F(p)
  loglik <- numeric(max_iter)
  criterion <- Inf
  count <- 0L
  monotone <- TRUE
  monotone_drops <- numeric(0)

  while (criterion > precision && count < max_iter) {
    count <- count + 1L
    if (verbose) {
      cat("Iteration:", count, "\n")
    }

    Sigma <- matrix_make_posdef(Sigma, eig_floor = eig_floor)
    Psi <- matrix_make_posdef(Psi, eig_floor = eig_floor)
    Sigma_inv <- solve(Sigma)
    Psi_inv <- solve(Psi)

    Wbar_list <- vector("list", n)
    S_list <- vector("list", n)
    wbar_vec_list <- vector("list", n)

    # E-step: W | X is N(Q^{-1}b, Q^{-1}) truncated to R_+^p.
    for (i in seq_len(n)) {
      Y <- X[, , i] - M
      b <- diag(A %*% Psi_inv %*% t(Y) %*% Sigma_inv)
      Q <- diag(p) + Sigma_inv * (A %*% Psi_inv %*% t(A))
      Q <- matrix_make_posdef(Q, eig_floor = eig_floor)
      Q_inv <- solve(Q)

      mt <- tmvtnorm::mtmvnorm(
        mean = as.vector(Q_inv %*% b),
        sigma = Q_inv,
        lower = rep(0, p),
        upper = rep(Inf, p)
      )

      wbar <- as.vector(mt$tmean)
      S <- as.matrix(mt$tvar)

      Wbar_list[[i]] <- diag(wbar, p, p)
      S_list[[i]] <- S
      wbar_vec_list[[i]] <- wbar
    }

    # CM-step for M: average of X_i - E(W_i | X_i) A.
    M_new <- matrix(0, p, q)
    for (i in seq_len(n)) {
      M_new <- M_new + X[, , i] - Wbar_list[[i]] %*% A
    }
    M_new <- M_new / n

    # CM-step for A.  S2 is E(w_i w_i' | X_i) = Var + mean mean'.
    Left <- matrix(0, p, p)
    Right <- matrix(0, p, q)
    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      wbar <- wbar_vec_list[[i]]
      S2 <- S_list[[i]] + wbar %*% t(wbar)
      Left <- Left + Sigma_inv * S2
      Right <- Right + Wbar_list[[i]] %*% Sigma_inv %*% Y
    }
    Left <- matrix_make_posdef(Left, eig_floor = eig_floor)
    A_new <- solve(Left, Right)

    # CM-step for Psi.  Delta is E(vec(Y - W A) vec(Y - W A)' | X).
    Psi_new <- matrix(0, q, q)
    G <- kronecker(t(A_new), diag(p)) %*% Fmat
    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      wbar <- wbar_vec_list[[i]]
      S2 <- S_list[[i]] + wbar %*% t(wbar)
      Delta <-
        matrix_vectorize(Y) %*% t(matrix_vectorize(Y)) -
        matrix_vectorize(Y) %*% t(matrix_vectorize(Wbar_list[[i]] %*% A_new)) -
        matrix_vectorize(Wbar_list[[i]] %*% A_new) %*% t(matrix_vectorize(Y)) +
        G %*% S2 %*% t(G)

      L <- mvrsn_chol_lower(matrix_make_posdef(Delta, eig_floor = eig_floor), "Delta")
      for (j in seq_len(p * q)) {
        Bij <- matrix(L[, j], nrow = p, ncol = q)
        Psi_new <- Psi_new + t(Bij) %*% Sigma_inv %*% Bij
      }
    }
    Psi_new <- matrix_make_posdef(Psi_new / (n * p), eig_floor = eig_floor)

    if (normalize_Psi) {
      det_Psi <- det(Psi_new)
      if (!is.finite(det_Psi) || det_Psi <= 0) {
        stop("Psi update has non-positive determinant.", call. = FALSE)
      }
      Psi_new <- matrix_make_posdef(Psi_new / det_Psi^(1 / q), eig_floor = eig_floor)
    }

    # CM-step for Sigma, now conditioning on the updated Psi.
    Psi_new_inv <- solve(Psi_new)
    Sigma_new <- matrix(0, p, p)
    for (i in seq_len(n)) {
      Y <- X[, , i] - M_new
      wbar <- wbar_vec_list[[i]]
      S2 <- S_list[[i]] + wbar %*% t(wbar)
      Delta <-
        matrix_vectorize(Y) %*% t(matrix_vectorize(Y)) -
        matrix_vectorize(Y) %*% t(matrix_vectorize(Wbar_list[[i]] %*% A_new)) -
        matrix_vectorize(Wbar_list[[i]] %*% A_new) %*% t(matrix_vectorize(Y)) +
        G %*% S2 %*% t(G)

      L <- mvrsn_chol_lower(matrix_make_posdef(Delta, eig_floor = eig_floor), "Delta")
      for (j in seq_len(p * q)) {
        Bij <- matrix(L[, j], nrow = p, ncol = q)
        Sigma_new <- Sigma_new + Bij %*% Psi_new_inv %*% t(Bij)
      }
    }
    Sigma_new <- matrix_make_posdef(Sigma_new / (n * q), eig_floor = eig_floor)

    M <- M_new
    A <- A_new
    Sigma <- Sigma_new
    Psi <- Psi_new

    loglik[count] <- loglik_mvrsn(X, M, A, Sigma, Psi)

    if (count > 1L) {
      diff <- loglik[count] - loglik[count - 1L]
      if (diff < -monotone_tol) {
        monotone <- FALSE
        monotone_drops <- c(monotone_drops, diff)
      }
      criterion <- abs(diff) / (abs(loglik[count - 1L]) + .Machine$double.eps)
    }
  }

  loglik <- loglik[seq_len(count)]

  if (count == max_iter && criterion > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - as.integer(normalize_Psi)
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(M = M,
                  A = A,
                  Sigma = Sigma,
                  Psi = Psi,
                  loglik = loglik[count],
                  loglik_history = loglik,
                  iterations = count,
                  BIC = BIC,
                  converged = (criterion <= precision),
                  criterion = criterion,
                  monotone = monotone,
                  monotone_drops = monotone_drops,
                  normalize_Psi = normalize_Psi,
                  npar = npar)

  class(obj.out) <- "MVRSN.ECM"
  obj.out

}


# ---------------------------------------------------------------------------
# Monte Carlo audit
# ---------------------------------------------------------------------------

#' Parameter set used in the paper's Monte Carlo section
#'
#' Constructs a predefined collection of `3` by `4` parameters for the MVRSN
#' Monte Carlo setting described in the row skew-normal article. The setting is
#' designed to study recovery of `M`, `A`, `Sigma`, and `Psi` when asymmetry is
#' driven by independent row-specific half-normal latent variables. `Psi` is
#' normalized to determinant one by default, matching the row/column scale
#' identifiability convention used by the implementation.
#'
#' @param normalize_Psi Logical scalar indicating whether `Psi` is normalized
#'   to determinant one, with reciprocal scaling applied to `Sigma`.
#' @return A named list containing `M`, `A`, `Sigma`, and `Psi`.
#' @details The parameter values are provided as a reproducible reference set;
#'   they are not fitted or estimated by this function.
#'
#' @examples
#' pars <- mvrsn_article_parameters()
#' names(pars)
#' det(pars$Psi)
#' @family MVCens MVRSN functions
#' @export
mvrsn_article_parameters <- function(normalize_Psi = TRUE) {
  M <- matrix(c(
    0.50, 1.00, 1.00, 1.00,
    1.50, 1.00, 1.50, 0.50,
    0.50, 1.00, 0.50, 1.50
  ), nrow = 3, ncol = 4)

  A <- matrix(c(
    0.18, 0.25, 0.71, 1.62,
    0.73, 1.42, 1.16, 0.83,
    0.22, 0.02, 0.13, 0.81
    ), nrow = 3, ncol = 4, byrow = TRUE)

  Sigma <- matrix(c(
    0.100, 0.040, 0.016,
    0.040, 0.100, 0.040,
    0.016, 0.040, 0.100
  ), nrow = 3, ncol = 3)

  Psi <- matrix(c(
    2.152, 1.721, 1.377, 1.102,
    1.721, 2.152, 1.721, 1.377,
    1.377, 1.721, 2.152, 1.721,
    1.102, 1.377, 1.721, 2.152
  ), nrow = 4, ncol = 4)

  if (normalize_Psi) {
    Psi <- Psi / det(Psi)^(1 / nrow(Psi))
  }

  list(M = M, A = A, Sigma = Sigma, Psi = Psi)
}

#' Summarize one ECM fit against known true parameters
#'
#' @param fit Return value from `mvrsn_ecm`.
#' @param truth List with true M, A, Sigma, and Psi.
#' @return A one-row data frame of errors and diagnostics.
#'
#' @keywords internal
#' @noRd
mvrsn_fit_diagnostics <- function(fit, truth) {
  data.frame(
    err_M = relative_frobenius_error(fit$M, truth$M),
    err_A = relative_frobenius_error(fit$A, truth$A),
    err_Sigma = relative_frobenius_error(fit$Sigma, truth$Sigma),
    err_Psi = relative_frobenius_error(fit$Psi, truth$Psi),
    loglik = fit$loglik,
    BIC = fit$BIC,
    iterations = fit$iterations,
    converged = fit$converged,
    monotone = fit$monotone,
    det_Psi = det(fit$Psi)
  )
}

#' Run a reproducible Monte Carlo audit of the MVRSN ECM estimator
#'
#' This function is the practical bridge between theory and evidence. It
#' repeatedly generates data from the stochastic representation, estimates the
#' parameters with the ECM algorithm, and measures whether estimates approach
#' the true values as n grows. For a quick smoke test use small values such as
#' `sample_sizes = c(10, 20)` and `replications = 2`; for paper-grade evidence
#' use the sample sizes and replication count described in the article.
#'
#' @param sample_sizes Integer vector of sample sizes.
#' @param replications Number of Monte Carlo replications per sample size.
#' @param truth Optional list with M, A, Sigma, and Psi. Defaults to article
#'   parameters.
#' @param max_iter Maximum ECM iterations for each fit.
#' @param precision ECM convergence tolerance.
#' @param seed Random seed.
#' @param verbose Print progress and ECM iterations.
#' @param workers Number of worker processes used for independent replications.
#'   Defaults to all logical CPU threads available to the host.
#' @return A list with per-replication results and a sample-size summary.
#'
#' @keywords internal
#' @noRd
mvrsn_monte_carlo <- function(sample_sizes = c(50, 100, 200, 400, 800, 1600),
                              replications = 200,
                              truth = mvrsn_article_parameters(),
                              max_iter = 200,
                              precision = 1e-6,
                              seed = 123,
                              verbose = FALSE,
                              workers = NULL) {
  if (any(sample_sizes <= 0L) || replications <= 0L) {
    stop("sample_sizes and replications must be positive.", call. = FALSE)
  }
  mvrsn_validate_parameters(M = truth$M, A = truth$A, Sigma = truth$Sigma, Psi = truth$Psi)

  tasks <- lapply(sample_sizes, function(n) {
    lapply(seq_len(replications), function(replication) {
      list(n = n, replication = replication)
    })
  })
  tasks <- unlist(tasks, recursive = FALSE)
  rows <- run_independent_tasks(tasks, function(task) {
      n <- task$n
      rep <- task$replication
      X <- rmvrsn(n, truth$M, truth$A, truth$Sigma, truth$Psi)
      fit <- tryCatch(
        mvrsn_ecm(X, max_iter = max_iter, precision = precision, verbose = FALSE),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        data.frame(
          n = n,
          replication = rep,
          err_M = NA_real_,
          err_A = NA_real_,
          err_Sigma = NA_real_,
          err_Psi = NA_real_,
          loglik = NA_real_,
          BIC = NA_real_,
          iterations = NA_integer_,
          converged = FALSE,
          monotone = FALSE,
          det_Psi = NA_real_,
          error = conditionMessage(fit)
        )
      } else {
        diag <- mvrsn_fit_diagnostics(fit, truth)
        cbind(
          data.frame(n = n, replication = rep),
          diag,
          data.frame(error = NA_character_)
        )
      }
  }, seed = seed, workers = workers)

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
      median_iterations = stats::median(df$iterations, na.rm = TRUE)
    )
  }))
  rownames(summary) <- NULL

  list(results = results, summary = summary, truth = truth)
}

mvcens_spec_mvrsn <- function() {
  new_model_spec(
    name = "MVRSN",
    validate = function(X = NULL, M = NULL, A = NULL, Sigma = NULL,
                        Psi = NULL, mode, ...) {
      if (mode == "fit") mvrsn_validate_sample_array(X)
      else mvrsn_validate_parameters(M = M, A = A, Sigma = Sigma, Psi = Psi)
      invisible(TRUE)
    },
    loglik = loglik_mvrsn,
    generate = function(n, M, A, Sigma, Psi, return_latent = FALSE, ...) {
      rmvrsn(n, M, A, Sigma, Psi, return_latent = return_latent)
    },
    parameter_count = function(p, q) model_parameter_count(p, q, skew = TRUE),
    fit = function(X, cc = NULL, LS = NULL, precision, max_iter,
                   normalize_Psi = TRUE, M_init = NULL, A_init = NULL,
                   Sigma_init = NULL, Psi_init = NULL, verbose = FALSE,
                   eig_floor = 1e-8, monotone_tol = 1e-7, ...) {
      mvrsn_ecm(X, max_iter = max_iter, precision = precision,
                normalize_Psi = normalize_Psi, M_init = M_init, A_init = A_init,
                Sigma_init = Sigma_init, Psi_init = Psi_init, verbose = verbose,
                eig_floor = eig_floor, monotone_tol = monotone_tol)
    }
  )
}
