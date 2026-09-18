#' Validate inputs for matrix-variate generators
#'
#' Checks whether the main inputs used by matrix-variate random generators have
#' compatible dimensions and valid types.
#'
#' @param n Positive integer. Number of matrix observations to generate.
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix.
#' @param Sigma Square numeric matrix of dimension \eqn{p \times p}. Row covariance matrix.
#' @param Psi Square numeric matrix of dimension \eqn{q \times q}. Column covariance matrix.
#' @param A Optional numeric matrix of dimension \eqn{p \times q}. Skewness matrix.
#' @param require_A Logical. If `TRUE`, the function requires `A` to be supplied.
#'
#' @return Invisibly returns `TRUE` if all inputs are valid.
#'
#' @keywords internal
validate_matrix_generator_inputs <- function(n, M, Sigma, Psi, A = NULL, require_A = FALSE) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.matrix(M)) {
    stop("'M' must be a matrix.")
  }

  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("'Sigma' must be a square matrix.")
  }

  if (!is.matrix(Psi) || nrow(Psi) != ncol(Psi)) {
    stop("'Psi' must be a square matrix.")
  }

  p <- nrow(Sigma)
  q <- ncol(Psi)

  if (!identical(dim(M), c(p, q))) {
    stop("'M' must have dimensions nrow(Sigma) x ncol(Psi).")
  }

  if (require_A) {
    if (is.null(A) || !is.matrix(A) || !identical(dim(A), c(p, q))) {
      stop("'A' must be a matrix with the same dimensions as 'M'.")
    }
  }

  invisible(TRUE)
}

rmvn_matrix <- function(n, M, Sigma, Psi) {
  validate_matrix_generator_inputs(n, M, Sigma, Psi)
  X <- array(NA_real_, dim = c(nrow(M), ncol(M), n))
  for (i in seq_len(n)) {
    X[, , i] <- LaplacesDemon::rmatrixnorm(M = M, U = Sigma, V = Psi)
  }
  X
}

#' Generate matrix-variate skew-normal samples
#'
#' Generates random matrices from a matrix-variate skew-normal-type model using
#' the stochastic representation
#'
#' \deqn{
#' X = A |Z| + E,
#' }
#'
#' where \eqn{Z \sim N(0,1)} and
#' \eqn{E \sim MN_{p \times q}(M, Sigma, Psi)}.
#'
#' @param n Positive integer. Number of matrices to generate.
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix used in the matrix-normal error.
#' @param A Numeric matrix of dimension \eqn{p \times q}. Skewness matrix.
#' @param Sigma Square numeric matrix of dimension \eqn{p \times p}. Row covariance matrix.
#' @param Psi Square numeric matrix of dimension \eqn{q \times q}. Column covariance matrix.
#'
#' @return A numeric array of dimension \eqn{p \times q \times n}, where each slice
#' `X[, , i]` is one simulated matrix observation.
#'
#' @keywords internal
rmvsn <- function(n, M, A, Sigma, Psi) {

  validate_matrix_generator_inputs(n = n, M = M, Sigma = Sigma, Psi = Psi, A = A, require_A = TRUE)

  p <- nrow(Sigma)
  q <- ncol(Psi)

  X <- array(NA_real_, dim = c(p, q, n))

  for (i in seq_len(n)) {
    z        <- abs(stats::rnorm(1))
    eps      <- LaplacesDemon::rmatrixnorm(M = M, U = Sigma, V = Psi)
    X[, , i] <- A * z + eps
  }

  X
}

#' Generate censored matrix-variate normal or skew-normal samples
#'
#' Simulates matrix-valued observations from either a matrix-variate normal or
#' skew-normal-type distribution and introduces censoring or missingness according
#' to a quantile-based threshold.
#'
#' The function first generates complete observations and then censors values below
#' the empirical `cens` quantile.
#'
#' @param n Positive integer. Number of matrices to generate.
#' @param cens Numeric scalar in `[0, 1]`. Proportion used to define the
#' censoring threshold. Zero returns complete data. One is supported for
#' interval censoring and mixed interval/missingness, but not for a wholly
#' missing sample (`Ind = 2`).
#' @param Ind Integer. Type of incomplete-data mechanism:
#' \itemize{
#'   \item `1`: interval censoring;
#'   \item `2`: missing values represented by `(-Inf, Inf)`;
#'   \item `3`: mixture of interval censoring and missing values.
#' }
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix.
#' @param Sigma Square numeric matrix of dimension \eqn{p \times p}. Row covariance matrix.
#' @param Psi Square numeric matrix of dimension \eqn{q \times q}. Column covariance matrix.
#' @param A Optional numeric matrix of dimension \eqn{p \times q}. Required when `dist = "SN"`.
#' @param dist Character. Distribution used to generate the complete data.
#' Must be either `"SN"` or `"Normal"`.
#'
#' @return A list with:
#' \describe{
#'   \item{X.cens}{Array containing the lower censoring bounds or observed values.}
#'   \item{cc}{Censoring indicator array. Entries equal to `1` indicate censored/missing values.}
#'   \item{LS}{Array containing upper bounds. Observed entries are set to zero.}
#' }
#'
#' @keywords internal
rmatrix_censored <- function(n, cens, Ind = 1, M, Sigma, Psi, A = NULL, dist = c("SN", "Normal")) {

  dist <- match.arg(dist)

  validate_matrix_generator_inputs(
    n = n,
    M = M,
    Sigma = Sigma,
    Psi = Psi,
    A = A,
    require_A = identical(dist, "SN")
  )

  if (!is.numeric(cens) || length(cens) != 1L || !is.finite(cens) || cens < 0 || cens > 1) {
    stop("'cens' must be a single number in the closed interval [0, 1].")
  }

  if (!is.numeric(Ind) || length(Ind) != 1L || !(Ind %in% c(1, 2, 3))) {
    stop("'Ind' must be one of 1, 2 or 3.")
  }

  if (cens == 1 && Ind == 2) {
    stop("'cens = 1' with 'Ind = 2' would create a wholly uninformative sample.")
  }

  p <- nrow(Sigma)
  q <- ncol(Psi)

  X.or <- array(NA_real_, dim = c(p, q, n))

  for (i in seq_len(n)) {
    eps <- LaplacesDemon::rmatrixnorm(M = M, U = Sigma, V = Psi)

    if (dist == "Normal") {
      X.or[, , i] <- eps
    } else {
      z <- abs(stats::rnorm(1))
      X.or[, , i] <- A * z + eps
    }
  }

  X.cens <- X.or
  LS <- X.or

  if (cens == 0) {
    cc <- array(0L, dim = dim(X.or))
  } else if (cens == 1) {
    cc <- array(1L, dim = dim(X.or))
  } else {
    cutoff <- as.numeric(stats::quantile(X.or, probs = cens, names = FALSE, type = 7))
    cc <- array(as.integer(X.or < cutoff), dim = dim(X.or))
  }

  cens_idx <- which(cc == 1L)
  obs_idx <- which(cc == 0L)

  LS[obs_idx] <- 0

  if (length(cens_idx) == 0L) {
    return(list(X.cens = X.cens, cc = cc, LS = LS))
  }

  if (Ind == 1) {
    spread <- stats::sd(X.or[cens_idx])

    if (!is.finite(spread) || spread <= 0) {
      spread <- .Machine$double.eps^0.5
    }

    LS[cens_idx] <- X.cens[cens_idx] + 2 * spread
    X.cens[cens_idx] <- X.cens[cens_idx] - 2 * spread
  }

  if (Ind == 2) {
    LS[cens_idx] <- Inf
    X.cens[cens_idx] <- -Inf
  }

  if (Ind == 3) {
    n_cens <- length(cens_idx)
    n_missing <- floor(n_cens / 2)
    n_interval <- n_cens - n_missing

    miss_idx <- cens_idx[seq_len(n_missing)]

    if (n_interval > 0L) {
      int_idx <- cens_idx[(n_missing + 1L):n_cens]
      spread <- stats::sd(X.or[int_idx])

      if (!is.finite(spread) || spread <= 0) {
        spread <- .Machine$double.eps^0.5
      }

      LS[int_idx] <- X.cens[int_idx] + 2 * spread
      X.cens[int_idx] <- X.cens[int_idx] - 2 * spread
    }

    if (n_missing > 0L) {
      LS[miss_idx] <- Inf
      X.cens[miss_idx] <- -Inf
    }
  }

  list(X.cens = X.cens, cc = cc, LS = LS)
}

#' Generate matrix-variate skew-t samples
#'
#' Generates random matrices from a matrix-variate skew-t-type distribution using
#' the stochastic representation
#'
#' \deqn{
#' X = M + \frac{A |Z| + E}{\sqrt{W}},
#' }
#'
#' where \eqn{Z \sim N(0,1)},
#' \eqn{E \sim MN_{p \times q}(0, Sigma, Psi)}, and
#' \eqn{W \sim Gamma(\nu/2, \nu/2)}.
#'
#' @param n Positive integer. Number of matrices to generate.
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix.
#' @param A Numeric matrix of dimension \eqn{p \times q}. Skewness matrix.
#' @param Sigma Square numeric matrix of dimension \eqn{p \times p}. Row covariance matrix.
#' @param Psi Square numeric matrix of dimension \eqn{q \times q}. Column covariance matrix.
#' @param nu Positive numeric scalar. Degrees-of-freedom parameter.
#'
#' @return A numeric array of dimension \eqn{p \times q \times n}.
#'
#' @keywords internal
rmvst <- function(n, M, A, Sigma, Psi, nu) {
  validate_matrix_generator_inputs(
    n = n,
    M = M,
    Sigma = Sigma,
    Psi = Psi,
    A = A,
    require_A = TRUE
  )

  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 0) {
    stop("'nu' must be a positive finite scalar.")
  }

  p <- nrow(Sigma)
  q <- ncol(Psi)

  X.or <- array(NA_real_, dim = c(p, q, n))
  M0 <- matrix(0, p, q)

  for (i in seq_len(n)) {
    z <- abs(stats::rnorm(1))
    w <- stats::rgamma(1, shape = nu / 2, rate = nu / 2)
    eps <- LaplacesDemon::rmatrixnorm(M = M0, U = Sigma, V = Psi)

    X.or[, , i] <- M + (A * z + eps) / sqrt(w)
  }

  X.or
}

#' Generate censored matrix-variate normal-inverse Gaussian samples
#'
#' Generates random matrices from a matrix-variate normal-inverse Gaussian-type
#' model using the stochastic representation
#'
#' \deqn{
#' X = M + W A + \sqrt{W} Z,
#' }
#'
#' where \eqn{Z \sim MN_{p \times q}(0, Sigma, Psi)} and
#' \eqn{W \sim IG(1, \tilde{\gamma})}.
#'
#' After generating the complete data, the function introduces censoring or
#' missingness based on the empirical `cens` quantile.
#'
#' @param n Positive integer. Number of matrices to generate.
#' @param cens Numeric scalar in `[0, 1]`. Proportion used to define the censoring threshold.
#' If `cens = 0`, no censoring is applied.
#' @param Ind Integer. Type of incomplete-data mechanism:
#' \itemize{
#'   \item `1`: interval censoring;
#'   \item `2`: missing values represented by `(-Inf, Inf)`;
#'   \item `3`: mixture of interval censoring and missing values.
#' }
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix.
#' @param Sigma Square numeric matrix of dimension \eqn{p \times p}. Row covariance matrix.
#' @param Psi Square numeric matrix of dimension \eqn{q \times q}. Column covariance matrix.
#' @param A Numeric matrix of dimension \eqn{p \times q}. Skewness matrix.
#' @param gamma_tilde Positive numeric scalar. Shape parameter of the inverse Gaussian latent variable.
#' Default is `2`.
#'
#' @return A list with:
#' \describe{
#'   \item{X.cens}{Array containing censored or observed data.}
#'   \item{cc}{Censoring indicator array. Entries equal to `1` indicate censored/missing values.}
#'   \item{LS}{Array containing upper censoring bounds.}
#'   \item{X.or}{Original complete generated data before censoring.}
#' }
#'
#' @keywords internal
rmvnig <- function(n, cens, Ind, M, Sigma, Psi, A, gamma_tilde = 2) {
  # ------------------------------------------------------------
  # Generate n samples from MVNIG_{p x q}(M, A, Sigma, Psi, gamma_tilde)
  #   X = M + W*A + sqrt(W)*Z,
  #   Z ~ MN_{p x q}(0, Sigma, Psi),  W ~ IG(mean = 1, shape = gamma_tilde)
  #
  # Ind:
  #   1 = interval censoring
  #   2 = missing
  #   3 = both
  # ------------------------------------------------------------

  if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
    stop("Package 'LaplacesDemon' is required.")
  }
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required.")
  }

  if (!is.numeric(n) || length(n) != 1L || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.numeric(cens) || length(cens) != 1L || cens < 0 || cens > 1) {
    stop("'cens' must be a proportion in [0, 1].")
  }

  if (!Ind %in% c(1, 2, 3)) {
    stop("'Ind' must be 1 (censoring), 2 (missing), or 3 (both).")
  }

  if (!is.matrix(M) || !is.matrix(Sigma) || !is.matrix(Psi) || !is.matrix(A)) {
    stop("'M', 'Sigma', 'Psi', and 'A' must be matrices.")
  }

  p <- nrow(Sigma)
  q <- nrow(Psi)

  if (ncol(Sigma) != p) {
    stop("'Sigma' must be a square matrix.")
  }
  if (ncol(Psi) != q) {
    stop("'Psi' must be a square matrix.")
  }
  if (!all(dim(M) == c(p, q))) {
    stop("'M' must have dimension p x q, where p = nrow(Sigma) and q = nrow(Psi).")
  }
  if (!all(dim(A) == c(p, q))) {
    stop("'A' must have dimension p x q, where p = nrow(Sigma) and q = nrow(Psi).")
  }

  if (!is.numeric(gamma_tilde) || length(gamma_tilde) != 1L || gamma_tilde <= 0) {
    stop("'gamma_tilde' must be a positive scalar.")
  }

  X.or <- array(NA_real_, dim = c(p, q, n))

  W <- statmod::rinvgauss(n, mean = 1, shape = gamma_tilde)

  Z0 <- matrix(0, p, q)
  for (i in seq_len(n)) {
    Z <- LaplacesDemon::rmatrixnorm(M = Z0, U = Sigma, V = Psi)
    X.or[, , i] <- M + W[i] * A + sqrt(W[i]) * Z
  }

  X.cens <- X.or
  LS <- array(0, dim = dim(X.or))
  cc <- array(0, dim = dim(X.or))

  if (cens == 0) {
    return(list(X.cens = X.cens, cc = cc, LS = LS, X.or = X.or))
  }

  cutoff <- as.numeric(stats::quantile(X.or, probs = cens, na.rm = TRUE, type = 7))
  cc <- (X.or < cutoff) + 0

  idx <- which(cc == 1)

  if (length(idx) == 0L) {
    return(list(X.cens = X.cens, cc = cc, LS = LS, X.or = X.or))
  }

  if (Ind == 1) {
    s_all <- stats::sd(X.or[idx])
    if (!is.finite(s_all) || s_all <= 0) {
      s_all <- 1
    }

    LS[idx] <- X.or[idx] + 2 * s_all
    X.cens[idx] <- X.or[idx] - 2 * s_all
  }

  if (Ind == 2) {
    LS[idx] <- Inf
    X.cens[idx] <- -Inf
  }

  if (Ind == 3) {
    n_idx <- length(idx)
    n_miss <- floor(0.5 * n_idx)
    miss_idx <- sample(idx, size = n_miss, replace = FALSE)
    cens_idx <- setdiff(idx, miss_idx)

    if (length(miss_idx) > 0L) {
      LS[miss_idx] <- Inf
      X.cens[miss_idx] <- -Inf
    }

    if (length(cens_idx) > 0L) {
      s_cens <- stats::sd(X.or[cens_idx])
      if (!is.finite(s_cens) || s_cens <= 0) {
        s_cens <- 1
      }

      LS[cens_idx] <- X.or[cens_idx] + 2 * s_cens
      X.cens[cens_idx] <- X.or[cens_idx] - 2 * s_cens
    }
  }

  list(
    X.cens = X.cens,
    cc = cc,
    LS = LS,
    X.or = X.or
  )
}

#' Generate random matrices from a matrix-variate variance-gamma distribution
#'
#' Generates random matrices from the stochastic representation
#' \deqn{
#'   X = M + W A + \sqrt{W} Psi,
#' }
#' where \eqn{W \sim \mathrm{Exp}(\mathrm{rate})} and
#' \eqn{Psi \sim \mathcal{N}_{p \times q}(0, \Sigma, \Psi)}.
#'
#' @param n Positive integer. Number of random matrices to generate.
#' @param M Numeric matrix of dimension \eqn{p \times q}. Location matrix.
#' @param A Numeric matrix of dimension \eqn{p \times q}. Skewness or drift matrix.
#' @param Sigma Square positive definite matrix. Row covariance matrix.
#' @param Psi Square positive definite matrix. Column covariance matrix.
#' @param rate Positive finite scalar. Rate parameter of the exponential
#' distribution used for \eqn{W}. Default is \code{1}.
#'
#' @return A numeric array of dimension
#' \eqn{p \times q \times number_of_samples}, where each slice
#' \code{X_array[, , i]} is one generated matrix.
#'
#' @details
#' For each sample, the function generates
#' \eqn{W \sim \mathrm{Exp}(\mathrm{rate})} and a standard normal matrix
#' \eqn{Z \in \mathbb{R}^{p \times q}}. Then it constructs
#' \eqn{Psi = L_\Sigma Z L_\Psi^\top}, where \eqn{L_\Sigma} and
#' \eqn{L_\Psi} are Cholesky factors of \eqn{\Sigma} and \eqn{\Psi},
#' respectively.
#'
#' The generated observation is then computed as
#' \deqn{
#'   X_i = M + W_i A + \sqrt{W_i} V_i.
#' }
#'
#' @keywords internal
rmvvg <- function(n, M, A, Sigma, Psi, rate = 1) {
  if (!is.numeric(n) ||
      length(n) != 1L ||
      !is.finite(n) ||
      n <= 0 ||
      n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.matrix(M)) {
    stop("'M' must be a matrix.")
  }

  if (!is.matrix(A)) {
    stop("'A' must be a matrix.")
  }

  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("'Sigma' must be a square matrix.")
  }

  if (!is.matrix(Psi) || nrow(Psi) != ncol(Psi)) {
    stop("'Psi' must be a square matrix.")
  }

  p <- nrow(M)
  q <- ncol(M)

  if (!all(dim(A) == c(p, q))) {
    stop("'A' must have the same dimensions as 'M'.")
  }

  if (!all(dim(Sigma) == c(p, p))) {
    stop("'Sigma' must have dimensions p x p, where p = nrow(M).")
  }

  if (!all(dim(Psi) == c(q, q))) {
    stop("'Psi' must have dimensions q x q, where q = ncol(M).")
  }

  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) || rate <= 0) {
    stop("'rate' must be a positive finite scalar.")
  }

  Sigma <- (Sigma + t(Sigma)) / 2
  Psi <- (Psi + t(Psi)) / 2

  chol_U <- tryCatch(
    chol(Sigma),
    error = function(e) stop("'Sigma' must be positive definite.")
  )

  chol_V <- tryCatch(
    chol(Psi),
    error = function(e) stop("'Psi' must be positive definite.")
  )

  X_array <- array(NA_real_, dim = c(p, q, n))

  for (i in seq_len(n)) {
    W <- stats::rexp(1, rate = rate)
    Z <- matrix(stats::rnorm(p * q), nrow = p, ncol = q)
    E <- chol_U %*% Z %*% t(chol_V)

    X_array[, , i] <- M + W * A + sqrt(W) * E
  }

  X_array
}

#' Generate random matrices from the MVRSN distribution
#'
#' This function implements the stochastic representation directly:
#' X_i = M + diag(w_i) A + V_i, where w_i has independent standard half-normal
#' entries and V_i is matrix-normal with covariance Psi %x% Sigma after
#' vectorization. It is the reference generator used by the Monte Carlo audit.
#'
#' @param n Number of matrices to generate.
#' @param M Location matrix of dimension p x q.
#' @param A Skewness matrix of dimension p x q.
#' @param Sigma Positive definite row covariance matrix.
#' @param Psi Positive definite column covariance matrix.
#' @param seed Optional random seed.
#' @param return_latent If TRUE, return the latent half-normal vectors too.
#' @return If `return_latent = FALSE`, an array p x q x n. Otherwise, a list
#'   with `X` and `W`.
#'
#' @keywords internal
rmvrsn <- function(n, M, A, Sigma, Psi, seed = NULL, return_latent = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(n) != 1L || n <= 0L) {
    stop("n must be a positive integer.", call. = FALSE)
  }

  pars <- mvrsn_validate_parameters(M = M, A = A, Sigma = Sigma, Psi = Psi)
  L_sigma <- mvrsn_chol_lower(pars$Sigma, "Sigma")
  R_psi <- chol(pars$Psi)

  X <- array(NA_real_, dim = c(pars$p, pars$q, n))
  W <- matrix(NA_real_, nrow = n, ncol = pars$p)

  for (i in seq_len(n)) {
    w <- abs(stats::rnorm(pars$p))
    Z <- matrix(stats::rnorm(pars$p * pars$q), nrow = pars$p, ncol = pars$q)
    noise <- L_sigma %*% Z %*% R_psi
    X[, , i] <- pars$M + diag(w, pars$p, pars$p) %*% pars$A  + noise
    W[i, ] <- w
  }

  if (return_latent) {
    return(list(X = X, W = W))
  }
  X
}

#' Generate random matrices from the MVREN distribution
#'
#' This function implements the stochastic representation directly:
#' \deqn{
#' X_i = M + \operatorname{diag}(w_i)A + V_i,
#' }
#' where the components of \eqn{w_i} are mutually independent exponential
#' random variables with rates given by `lambda`, and \eqn{V_i} follows a
#' matrix-normal distribution with covariance `Psi \%x\% Sigma` after
#' vectorization. It is the reference generator used by the Monte Carlo audit.
#'
#' @param n Number of matrices to generate.
#' @param M Location matrix of dimension `p` by `q`.
#' @param A Row-specific latent-effect matrix of dimension `p` by `q`.
#' @param Sigma Positive-definite row covariance matrix of dimension `p` by
#'   `p`.
#' @param Psi Positive-definite column covariance matrix of dimension `q` by
#'   `q`.
#' @param lambda Optional positive rate vector of length `p`. When `NULL`, a
#'   vector of ones is used. The high-level `mv_random("MVREN", ...)` API
#'   accepts only the unit vector; non-unit values here are retained for the
#'   unexported general stochastic representation.
#' @param return_latent If `TRUE`, return the generated latent exponential
#'   vectors together with the observations.
#'
#' @return If `return_latent = FALSE`, an array of dimension `p` by `q` by `n`.
#'   Otherwise, a list containing `X` and `W`, where `W` is an `n` by `p`
#'   matrix of latent exponential variables.
#'
#' @keywords internal
rmvren <- function(n,
                   M,
                   A,
                   Sigma,
                   Psi,
                   lambda = NULL,
                   return_latent = FALSE) {
  if (length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  pars <- mvren_validate_parameters(
    M = M,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    lambda = lambda
  )

  L_Sigma <- mvren_chol_lower(pars$Sigma, "Sigma")
  L_Psi <- mvren_chol_lower(pars$Psi, "Psi")

  X <- array(0, dim = c(pars$p, pars$q, n))
  W <- matrix(0, nrow = n, ncol = pars$p)

  for (i in seq_len(n)) {
    w <- stats::rexp(pars$p, rate = pars$lambda)
    Z <- matrix(stats::rnorm(pars$p * pars$q), nrow = pars$p, ncol = pars$q)
    noise <- L_Sigma %*% Z %*% t(L_Psi)

    W[i, ] <- w
    X[, , i] <- pars$M + diag(w, pars$p, pars$p) %*% pars$A  + noise
  }

  if (isTRUE(return_latent)) {
    return(list(X = X, W = W))
  }
  X
}
