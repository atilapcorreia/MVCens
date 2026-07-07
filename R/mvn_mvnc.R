#' Make a symmetric matrix positive definite
#'
#' Symmetrizes a square matrix and, if necessary, shifts its eigenvalues so that
#' the smallest eigenvalue is greater than `epsilon`.
#'
#' @param X Numeric square matrix.
#' @param epsilon Positive tolerance used as the minimum eigenvalue threshold.
#'
#' @return A symmetric positive definite matrix.
#'
#' @keywords internal
make_posdef <- function(X, epsilon = 1e-8) {

  X <- (X + t(X)) / 2

  eig <- eigen(X, symmetric = TRUE)

  lambda_min <- min(eig$values)

  if (lambda_min > epsilon) {

    return(X)

  }

  shift <- epsilon - lambda_min

  X + shift * diag(nrow(X))
}

#' Vectorize a matrix by columns
#'
#' Converts a matrix into a vector using column-wise ordering.
#'
#' @param X Numeric matrix.
#'
#' @return A numeric vector containing the entries of `X` stacked by columns.
#'
#' @keywords internal
vec_col <- function(X) {

  as.vector(matrixNormal::vec(X))

}

#' Invert a symmetric positive definite matrix
#'
#' Stabilizes a matrix to be positive definite and computes its inverse using
#' the Cholesky decomposition.
#'
#' @param A Numeric square matrix.
#' @param epsilon Positive tolerance used for positive-definiteness correction.
#'
#' @return The inverse of `A`.
#'
#' @keywords internal
solve_sym_pd <- function(A, epsilon = 1e-8) {

  A <- make_posdef((A + t(A)) / 2, epsilon = epsilon)

  chol2inv(chol(A))

}

#' Normalize a covariance matrix under determinant constraint
#'
#' Normalizes a covariance matrix so that its determinant is equal to one. This
#' is useful for identifiability in matrix-variate models, where the covariance
#' structure is represented through row and column scale matrices.
#'
#' @param A Numeric square covariance matrix.
#' @param epsilon Positive tolerance used for positive-definiteness correction.
#'
#' @return A positive definite matrix with determinant approximately equal to one.
#'
#' @keywords internal
normalize_cov_constraint <- function(A, epsilon = 1e-8) {

  A <- make_posdef((A + t(A)) / 2, epsilon = epsilon)

  logdetA <- as.numeric(determinant(A, logarithm = TRUE)$modulus)

  A / exp(logdetA / ncol(A))

}

#' Compute the ECM convergence criterion
#'
#' Computes an Aitken-type convergence criterion based on the last three
#' log-likelihood values.
#'
#' @param loglik Numeric vector containing the log-likelihood values.
#' @param iter Current iteration index.
#'
#' @return A nonnegative scalar convergence criterion.
#'
#' @keywords internal
compute_ecm_criterion <- function(loglik, iter) {

  if (iter <= 1L) {
    return(Inf)
  }

  ll_curr <- loglik[iter]
  ll_prev <- loglik[iter - 1L]

  if (!is.finite(ll_curr) || !is.finite(ll_prev)) {
    return(Inf)
  }

  abs(1 - ll_prev / ll_curr)
}

#' Validate matrix-variate normal input data
#'
#' Checks whether the sample object is a three-dimensional array without missing
#' values.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#'
#' @return Invisibly returns `TRUE` if the input is valid.
#'
#' @keywords internal
validate_mvn_input <- function(samples) {

  if (length(dim(samples)) != 3L) {

    stop("'samples' deve ser um array 3D com dimensoes p x q x n.")

  }

  if (anyNA(samples)) {

    stop("'samples' contem NA. Trate os dados antes de rodar o ECM.")

  }

  invisible(TRUE)

}

#' Validate censored matrix-variate normal input data
#'
#' Checks whether the observed/censored samples, censoring indicators, and upper
#' limits have compatible dimensions and valid entries.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param cc Integer array with the same dimensions as `samples`. Entries equal
#' to `1` indicate censored or missing observations, and entries equal to `0`
#' indicate observed values.
#' @param LS Numeric array with the same dimensions as `samples`, containing
#' upper censoring limits.
#'
#' @return Invisibly returns `TRUE` if all inputs are valid.
#'
#' @keywords internal
validate_mvnc_input <- function(samples, cc, LS) {

  validate_mvn_input(samples)

  if (!identical(dim(samples), dim(cc)) || !identical(dim(samples), dim(LS))) {

    stop("'samples', 'cc' e 'LS' devem ter as mesmas dimensoes.")

  }

  if (anyNA(cc) || anyNA(LS)) {

    stop("'cc' e 'LS' nao podem conter NA.")

  }

  if (!all(cc %in% c(0, 1))) {

    stop("'cc' deve conter apenas 0 e 1.")

  }

  invisible(TRUE)

}

#' Compute a fill value for censored or missing entries
#'
#' Computes a stable initial value used to replace censored or missing entries
#' before starting the ECM algorithm.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param cc Optional censoring indicator array. If supplied, entries with
#' `cc == 0` are treated as observed.
#'
#' @return A finite numeric scalar.
#'
#' @keywords internal
get_observed_fill_value <- function(samples, cc = NULL) {

  observed <- if (is.null(cc)) samples else samples[cc == 0]

  fill_value <- mean(observed)

  if (!is.finite(fill_value)) {

    fill_value <- mean(samples)

  }

  if (!is.finite(fill_value)) {

    fill_value <- 0

  }

  fill_value

}

#' Compute moments of a truncated multivariate normal vector
#'
#' Computes the first and second conditional moments required in the E-step of
#' the ECM algorithm for censored matrix-variate normal data.
#'
#' @param y Numeric vector representing one vectorized matrix observation.
#' @param cc1 Integer vector indicating censored entries. Entries equal to `1`
#' are treated as censored/missing.
#' @param LS1 Numeric vector of upper censoring limits.
#' @param mu1 Mean vector.
#' @param Vari Covariance matrix of the vectorized matrix observation.
#' @param epsilon Positive tolerance used for numerical stabilization.
#'
#' @return A list with components:
#' \describe{
#'   \item{tuy}{Conditional first moment \eqn{E(Y \mid data)}.}
#'   \item{tuyy}{Conditional second moment \eqn{E(YY^\top \mid data)}.}
#'   \item{censored}{Logical value indicating whether censoring was present.}
#' }
#'
#' @keywords internal
get_truncated_moments <- function(y, cc1, LS1, mu1, Vari, epsilon = 1e-8) {

  miss <- which(cc1 == 1)

  obs <- which(cc1 == 0)

  d <- length(y)

  if (length(miss) == 0L) {

    tuy <- matrix(y, d, 1)

    return(list(tuy = tuy, tuyy = tcrossprod(tuy), censored = FALSE))

  }

  if (length(miss) == d) {

    aux <- MomTrunc::meanvarTMD(y, LS1, mu1, Vari, dist = "normal")

    tuy <- matrix(aux$mean, d, 1)

    tuyy <- if (!is.null(aux$EYY)) aux$EYY else aux$varcov + tcrossprod(tuy)

    return(list(tuy = tuy, tuyy = tuyy, censored = TRUE))

  }

  V_oo <- Vari[obs, obs, drop = FALSE]
  V_mo <- Vari[miss, obs, drop = FALSE]
  V_om <- Vari[obs, miss, drop = FALSE]
  V_mm <- Vari[miss, miss, drop = FALSE]

  V_oo_inv <- solve_sym_pd(V_oo, epsilon = epsilon)
  muc <- mu1[miss] + V_mo %*% V_oo_inv %*% (y[obs] - mu1[obs])
  Sc <- V_mm - V_mo %*% V_oo_inv %*% V_om
  Sc <- make_posdef((Sc + t(Sc)) / 2, epsilon = epsilon)

  aux <- MomTrunc::meanvarTMD(y[miss], LS1[miss], as.vector(muc), Sc, dist = "normal")

  tuy <- matrix(y, d, 1)
  tuy[miss] <- aux$mean

  tuyy <- matrix(0, d, d)
  tuyy[obs, obs] <- tcrossprod(y[obs])
  tuyy[obs, miss] <- tcrossprod(y[obs], as.vector(tuy[miss]))
  tuyy[miss, obs] <- tcrossprod(as.vector(tuy[miss]), y[obs])
  tuyy[miss, miss] <- aux$varcov + tcrossprod(as.vector(tuy[miss]))

  list(tuy = tuy, tuyy = tuyy, censored = TRUE)

}

#' Convert truncated moments into a covariance-type matrix
#'
#' Computes the conditional second central moment matrix from the first and
#' second truncated moments.
#'
#' @param tuy Conditional first moment vector.
#' @param tuyy Conditional second moment matrix.
#' @param mu1 Current mean vector.
#' @param epsilon Positive tolerance used for numerical stabilization.
#'
#' @return A positive definite covariance-type matrix.
#'
#' @keywords internal
moment_to_omega <- function(tuy, tuyy, mu1, epsilon = 1e-8) {

  omega <- tuyy - tuy %*% t(mu1) - mu1 %*% t(tuy) + mu1 %*% t(mu1)
  omega <- as.matrix(Matrix::nearPD((omega + t(omega)) / 2)$mat)
  make_posdef(omega, epsilon = epsilon)

}

#' Compute censored sufficient statistics for one matrix observation
#'
#' Vectorizes one matrix observation, computes truncated moments when needed,
#' and returns the corresponding conditional statistics used in the ECM updates.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param cc Censoring indicator array.
#' @param LS Array of upper censoring limits.
#' @param mu Current estimate of the location matrix.
#' @param Vari Covariance matrix of the vectorized observation.
#' @param j Integer index of the matrix observation.
#' @param epsilon Positive tolerance used for numerical stabilization.
#'
#' @return A list containing vectorized data, censoring indicators, conditional
#' moments, predicted mean matrix, and covariance-type correction matrix.
#'
#' @keywords internal
get_censored_sample_stats <- function(samples, cc, LS, mu, Vari, j, epsilon = 1e-8) {

  y1 <- vec_col(samples[, , j])
  cc1 <- vec_col(cc[, , j])
  LS1 <- vec_col(LS[, , j])
  mu1 <- matrix(vec_col(mu), ncol = 1)

  moments <- get_truncated_moments(
    y = y1,
    cc1 = cc1,
    LS1 = LS1,
    mu1 = as.vector(mu1),
    Vari = Vari,
    epsilon = epsilon
  )

  list(
    cc1 = cc1,
    y1 = y1,
    LS1 = LS1,
    mu1 = mu1,
    mean_mat = matrix(moments$tuy, nrow(mu), ncol(mu)),
    tuy = moments$tuy,
    tuyy = moments$tuyy,
    omega = if (moments$censored) moment_to_omega(moments$tuy, moments$tuyy, mu1, epsilon = epsilon) else NULL
  )
}

#' Matrix-variate normal log-likelihood
#'
#' Computes the observed log-likelihood for complete matrix-variate normal data.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param M Location matrix of dimension \eqn{p \times q}.
#' @param Sigma Row covariance matrix of dimension \eqn{p \times p}.
#' @param Psi Column covariance matrix of dimension \eqn{q \times q}.
#'
#' @return Numeric scalar containing the total log-likelihood.
#'
#' @export
loglik_mvn <- function(samples, M, Sigma, Psi) {

  validate_mvn_input(samples)

  n <- dim(samples)[3]

  auxiliary_variable <- 0

  for (j in seq_len(n)) {
    auxiliary_variable <- auxiliary_variable + matrixNormal::dmatnorm(samples[, , j],
                                                                      M,
                                                                      Sigma,
                                                                      Psi,
                                                                      tol = .Machine$double.eps^0.5,
                                                                      log = TRUE)
  }

  auxiliary_variable

}

#' Censored matrix-variate normal log-likelihood
#'
#' Computes the observed log-likelihood for censored matrix-variate normal data.
#' Fully observed matrices are evaluated through the matrix-normal density.
#' Fully censored matrices are evaluated through the matrix-normal probability,
#' and partially censored matrices are handled by conditioning on the observed
#' entries.
#'
#' @param cc Censoring indicator array. Entries equal to `1` indicate censored or
#' missing values.
#' @param LS Array of upper censoring limits.
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param M Location matrix.
#' @param Sigma Row covariance matrix.
#' @param Psi Column covariance matrix.
#' @param epsilon Positive tolerance used for numerical stabilization.
#'
#' @return Numeric scalar containing the censored-data log-likelihood.
#'
#' @export
loglik_mvnc <- function(cc, LS, samples, M, Sigma, Psi, epsilon = 1e-8) {

  validate_mvnc_input(samples, cc, LS)

  p <- dim(samples)[1]
  q <- dim(samples)[2]
  n <- dim(samples)[3]

  ver <- numeric(n)
  Vari <- kronecker(Psi, Sigma)
  mu1 <- vec_col(M)

  for (j in seq_len(n)) {

    cc1  <- vec_col(cc[, , j])
    LS1  <- vec_col(LS[, , j])
    y1   <- vec_col(samples[, , j])
    miss <- which(cc1 == 1)
    obs  <- which(cc1 == 0)

    if (length(miss) == 0L) {

      ver[j] <- matrixNormal::dmatnorm(samples[, , j],
                                      M,
                                      Sigma,
                                      Psi,
                                      tol = .Machine$double.eps^0.5,
                                      log = TRUE)
      next
    }

    if (length(miss) == p * q) {

      ver[j] <- matrixNormal::pmatnorm(Lower = samples[, , j],
                                       Upper = LS[, , j],
                                       M,
                                       Sigma,
                                       Psi,
                                       tol = .Machine$double.eps^0.5,
                                       log = TRUE)
      next
    }

    V_oo <- make_posdef((Vari[obs, obs, drop = FALSE] + t(Vari[obs, obs, drop = FALSE])) / 2, epsilon = epsilon)
    V_mo <- Vari[miss, obs, drop = FALSE]
    V_om <- Vari[obs, miss, drop = FALSE]
    V_mm <- Vari[miss, miss, drop = FALSE]
    V_oo_inv <- solve_sym_pd(V_oo, epsilon = epsilon)

    muc <- mu1[miss] + V_mo %*% V_oo_inv %*% (y1[obs] - mu1[obs])
    Sc <- V_mm - V_mo %*% V_oo_inv %*% V_om
    Sc <- make_posdef((Sc + t(Sc)) / 2, epsilon = epsilon)

    auxcdf <- mnormt::sadmvn(lower = as.vector(y1[miss]),
                             upper = as.vector(LS1[miss]),
                             mean = as.vector(muc),
                             varcov = Sc,
                             abseps = .Machine$double.eps^0.5)

    if (!is.finite(auxcdf) || auxcdf <= 0) {

      auxcdf <- .Machine$double.xmin

    }

    ver[j] <- mnormt::dmnorm(x = as.vector(y1[obs]),
                             mean = as.vector(mu1[obs]),
                             varcov = V_oo,
                             log = TRUE) + log(auxcdf)
  }

  sum(ver)

}

#' Compute covariance update components
#'
#' Computes the matrix sums required for updating the row and column covariance
#' matrices in the ECM algorithm.
#'
#' @param L1 Matrix whose columns represent vectorized residual-like quantities.
#' @param Sigma Current row covariance matrix.
#' @param Psi Current column covariance matrix.
#' @param epsilon Positive tolerance used for numerical stabilization.
#'
#' @return A list with:
#' \describe{
#'   \item{sPsi}{Accumulated component used to update `Psi`.}
#'   \item{sSigma}{Accumulated component used to update `Sigma`.}
#' }
#'
#' @keywords internal
somaL3 <- function(L1, Sigma, Psi, epsilon = 1e-8) {

  n <- ncol(L1)
  p <- ncol(Sigma)
  q <- ncol(Psi)

  Sigma_inv <- solve_sym_pd(Sigma, epsilon = epsilon)
  Psi_inv <- solve_sym_pd(Psi, epsilon = epsilon)

  suma2 <- matrix(0, q, q)
  suma3 <- matrix(0, p, p)

  for (j in seq_len(n)) {

    aux <- matrix(L1[, j], p, q)
    suma2 <- suma2 + t(aux) %*% Sigma_inv %*% aux
    suma3 <- suma3 + aux %*% Psi_inv %*% t(aux)

  }

  list(sPsi = suma2, sSigma = suma3)
}

#' ECM estimation for the censored matrix-variate normal model
#'
#' Fits a censored matrix-variate normal model by an ECM algorithm. The model
#' assumes matrix-valued observations
#'
#' \deqn{
#' X_i \sim MN_{p \times q}(\mu, \Sigma, \Psi),
#' }
#'
#' possibly subject to censoring or missingness. The algorithm alternates between
#' computing conditional moments of censored entries and updating the location,
#' row covariance, and column covariance matrices.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param cc Censoring indicator array with the same dimensions as `samples`.
#' Entries equal to `1` indicate censored or missing values.
#' @param LS Array of upper censoring limits.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#'
#' @return An object of class `"MVNC.ECM"` containing:
#' \describe{
#'   \item{mu}{Estimated location matrix.}
#'   \item{Sigma}{Estimated row covariance matrix.}
#'   \item{Psi}{Estimated column covariance matrix.}
#'   \item{dadosPred}{Array with predicted conditional values for censored entries.}
#'   \item{loglik}{Final observed log-likelihood.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iter}{Number of iterations performed.}
#'   \item{converged}{Logical value indicating whether convergence was reached.}
#' }
#'
#' @keywords internal
mvnc_ecm <- function(samples, cc, LS, precision = 1e-7, MaxIter = 50) {

  validate_mvnc_input(samples, cc, LS)

  p <- dim(samples)[1]
  q <- dim(samples)[2]
  n <- dim(samples)[3]

  loglik <- numeric(MaxIter)
  criterio <- Inf
  count <- 0L

  dadosPred <- samples
  dadosInit <- samples
  dadosInit[cc == 1] <- get_observed_fill_value(samples, cc)

  mu    <- apply(dadosInit, c(1, 2), mean)
  Psi   <- diag(q)
  Sigma <- diag(p)
  Vari  <- kronecker(Psi, Sigma)

  while (criterio > precision && count < MaxIter) {

    count <- count + 1L

    suma0 <- matrix(0, p, q)
    suma2 <- matrix(0, q, q)
    suma3 <- matrix(0, p, p)

    for (j in seq_len(n)) {
      stats_mu         <- get_censored_sample_stats(samples, cc, LS, mu, Vari, j)
      dadosPred[, , j] <- stats_mu$mean_mat
      suma0            <- suma0 + stats_mu$mean_mat
    }

    mu <- suma0 / n
    Sigma_inv <- solve_sym_pd(Sigma)

    for (j in seq_len(n)) {

      stats_psi <- get_censored_sample_stats(samples, cc, LS, mu, Vari, j)

      if (sum(stats_psi$cc1) == 0L) {

        Xc    <- stats_psi$mean_mat - mu
        suma2 <- suma2 + t(Xc) %*% Sigma_inv %*% Xc

      } else {

        aux   <- somaL3(t(chol(stats_psi$omega)), Sigma, Psi)
        suma2 <- suma2 + aux$sPsi
        dadosPred[, , j] <- stats_psi$mean_mat

      }
    }

    Psi     <- normalize_cov_constraint(suma2)
    Vari    <- kronecker(Psi, Sigma)
    Psi_inv <- solve_sym_pd(Psi)

    for (j in seq_len(n)) {

      stats_sigma <- get_censored_sample_stats(samples, cc, LS, mu, Vari, j)

      if (sum(stats_sigma$cc1) == 0L) {

        Xc    <- stats_sigma$mean_mat - mu
        suma3 <- suma3 + Xc %*% Psi_inv %*% t(Xc)

      } else {

        aux   <- somaL3(t(chol(stats_sigma$omega)), Sigma, Psi)
        suma3 <- suma3 + aux$sSigma
        dadosPred[, , j] <- stats_sigma$mean_mat

      }
    }

    Sigma <- make_posdef((suma3 / (q * n) + t(suma3 / (q * n))) / 2)
    Vari  <- kronecker(Psi, Sigma)
    loglik[count] <- loglik_mvnc(cc, LS, samples, mu, Sigma, Psi)
    criterio <- compute_ecm_criterion(loglik, count)
  }

  loglik <- loglik[seq_len(count)]

  if (count == MaxIter && criterio > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(mu = mu,
                  Sigma = Sigma,
                  Psi = Psi,
                  dadosPred = dadosPred,
                  loglik = loglik[count],
                  BIC = BIC,
                  iter = count,
                  converged = (criterio <= precision))

  class(obj.out) <- "MVNC.ECM"
  obj.out
}

#' ECM estimation for the matrix-variate normal model
#'
#' Fits a complete-data matrix-variate normal model by an ECM algorithm. The model
#' assumes
#'
#' \deqn{
#' X_i \sim MN_{p \times q}(\mu, \Sigma, \Psi),
#' }
#'
#' where `mu` is the location matrix, `Sigma` is the row covariance matrix, and
#' `Psi` is the column covariance matrix.
#'
#' @param samples Numeric array with dimensions \eqn{p \times q \times n}.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#'
#' @return An object of class `"MVN.ECM"` containing:
#' \describe{
#'   \item{mu}{Estimated location matrix.}
#'   \item{Sigma}{Estimated row covariance matrix.}
#'   \item{Psi}{Estimated column covariance matrix.}
#'   \item{loglik}{Final log-likelihood.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iter}{Number of iterations performed.}
#'   \item{converged}{Logical value indicating whether convergence was reached.}
#' }
#'
#' @keywords internal
mvn_ecm <- function(samples, precision = 1e-7, MaxIter = 50) {

  validate_mvn_input(samples)

  p <- dim(samples)[1]
  q <- dim(samples)[2]
  n <- dim(samples)[3]

  loglik   <- numeric(MaxIter)
  criterio <- Inf
  count    <- 0L

  mu    <- apply(samples, c(1, 2), mean)
  Psi   <- diag(q)
  Sigma <- diag(p)

  while (criterio > precision && count < MaxIter) {

    count <- count + 1L

    Sigma_inv <- solve_sym_pd(Sigma)
    suma2     <- matrix(0, q, q)

    for (j in seq_len(n)) {

      Xc <- samples[, , j] - mu
      suma2 <- suma2 + t(Xc) %*% Sigma_inv %*% Xc

    }

    Psi     <- normalize_cov_constraint(suma2)
    Psi_inv <- solve_sym_pd(Psi)
    suma3   <- matrix(0, p, p)

    for (j in seq_len(n)) {

      Xc <- samples[, , j] - mu
      suma3 <- suma3 + Xc %*% Psi_inv %*% t(Xc)

    }

    Sigma         <- make_posdef((suma3 / (q * n) + t(suma3 / (q * n))) / 2)
    loglik[count] <- loglik_mvn(samples, mu, Sigma, Psi)
    criterio      <- compute_ecm_criterion(loglik, count)

  }

  loglik <- loglik[seq_len(count)]

  if (count == MaxIter && criterio > precision) {

    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")

  }

  npar <- (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(mu = mu,
                  Sigma = Sigma,
                  Psi = Psi,
                  loglik = loglik[count],
                  BIC = BIC,
                  iter = count,
                  converged = (criterio <= precision))

  class(obj.out) <- "MVN.ECM"
  obj.out
}
