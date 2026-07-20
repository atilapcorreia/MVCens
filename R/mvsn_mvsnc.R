#' Matrix square root
#'
#' Computes the symmetric square root of a nonnegative scalar or symmetric
#' positive semidefinite matrix.
#'
#' @param A Numeric scalar or numeric square matrix.
#' @param epsilon Numerical tolerance for detecting negative eigenvalues.
#'
#' @return A matrix containing the square root of `A`.
#'
#' @keywords internal
#' @noRd
sqrtm <- function(A, epsilon = 1e-8) {

  if (length(A) == 1L) {

    if (!is.numeric(A) || !is.finite(A) || A < 0) {
      stop("Matrix square root is not defined.")

    }

    return(matrix(sqrt(A), 1, 1))

  }

  A <- as.matrix(A)

  if (!is.numeric(A) || nrow(A) != ncol(A)) {
    stop("'A' must be a numeric square matrix.")
  }

  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE)

  if (any(!is.finite(eig$values))) {
    stop("Matrix square root could not be computed.")
  }

  if (min(eig$values) < -epsilon) {
    stop("Matrix square root is not defined for matrices with negative eigenvalues.")
  }

  vals <- pmax(eig$values, 0)
  Asqrt <- eig$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(eig$vectors)
  as.matrix((Asqrt + t(Asqrt)) / 2)

  }

#' Inverse Mills ratio
#'
#' Computes the inverse Mills ratio for a normal random variable.
#'
#' @param x Finite numeric scalar.
#' @param mu Mean of the normal distribution. Default is `0`.
#' @param sd Standard deviation of the normal distribution. Default is `1`.
#'
#' @return A positive numeric scalar.
#'
#' @keywords internal
#' @noRd
invmills <- function(x, mu = 0, sd = 1) {

  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("'x' must be a finite scalar.")
  }

  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) {
    stop("'mu' must be a finite scalar.")
  }

  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("'sd' must be a positive finite scalar.")
  }

  z <- (x - mu) / sd

  if (z < -1e4) {
    return(-z / sd)
  }

  exp(stats::dnorm(x, mean = mu, sd = sd, log = TRUE) - stats::pnorm(q = x, mean = mu, sd = sd, log.p = TRUE))

}

#' Multivariate normal rectangle probability
#'
#' Computes the probability that a multivariate normal random vector lies inside
#' a rectangular region. For low-dimensional problems, `mvtnorm::pmvnorm()` is
#' tried first. Otherwise, or if it fails, `tlrmvnmvt::pmvn()` is used.
#'
#' @param lower Numeric vector of lower integration limits.
#' @param upper Numeric vector of upper integration limits.
#' @param mean Numeric mean vector.
#' @param sigma Positive definite covariance matrix.
#' @param nu Optional degrees of freedom. Currently only `nu = NULL` is supported.
#' @param uselog2 Logical. If `TRUE`, returns the base-2 logarithm of the probability.
#' @param epsilon Numerical tolerance for covariance regularization.
#'
#' @return A numeric probability or its base-2 logarithm.
#'
#' @keywords internal
#' @noRd
prob_opt <- function(lower = rep(-Inf, ncol(sigma)),
                     upper = rep(Inf, ncol(sigma)),
                     mean = rep(0, ncol(sigma)),
                     sigma,
                     nu = NULL,
                     uselog2 = FALSE,
                     epsilon = 1e-8) {

  sigma <- make_posdef(as.matrix(sigma), epsilon = epsilon)
  p <- ncol(sigma)
  mean <- c(mean)
  lower <- c(lower)
  upper <- c(upper)

  if (nrow(sigma) != p) {

    stop("'sigma' must be a square matrix.")

  }

  if (length(mean) != p || length(lower) != p || length(upper) != p) {

    stop("'mean', 'lower' and 'upper' must have length equal to ncol(sigma).")

  }

  if (!is.null(nu)) {

    stop("Only the normal case is implemented here ('nu = NULL').")

  }

  if (p < 10L) {
    prob <- tryCatch(mvtnorm::pmvnorm(lower = lower,
                                      upper = upper,
                                      mean = mean,
                                      sigma = sigma)[1],
                                      error = function(e) NA_real_
    )

    if (is.finite(prob) && prob >= 0) {
      return(if (uselog2) log2(prob) else prob)
    }
  }

  tlrmvnmvt::pmvn(lower = lower,
                  upper = upper,
                  mean = mean,
                  sigma = sigma,
                  uselog2 = uselog2)[[1]]
}

#' Density of the multivariate skew-normal distribution
#'
#' Evaluates the probability density of a multivariate skew-normal distribution
#' at a given observation, with numerical regularization of the covariance
#' matrices controlled by \code{epsilon}.
#'
#' @param y Numeric vector containing the values at which the density is
#'   evaluated.
#' @param mu Numeric vector specifying the location parameter. It must have the
#'   same length as \code{y}.
#' @param Sigma Numeric positive definite covariance matrix associated with the
#'   symmetric component of the distribution. Its number of rows and columns
#'   must equal the length of \code{mu}.
#' @param lambda Numeric vector specifying the direction and magnitude of
#'   skewness. It must have the same length as \code{mu}.
#' @param epsilon Positive numeric scalar specifying the tolerance used for
#'   covariance-matrix regularization and numerical positive-definiteness
#'   checks. The default is \code{1e-8}.
#'
#' @return A numeric scalar containing the density evaluated at \code{y}.
#'
#' @export

dmvsn <- function(y, mu, Sigma, lambda, epsilon = 1e-8) {

  y      <- c(y)
  mu     <- c(mu)
  Sigma  <- make_posdef(as.matrix(Sigma), epsilon = epsilon)
  lambda <- matrix(c(lambda), ncol = 1)

  p <- length(mu)

  if (length(y) != p) {
    stop("'y' and 'mu' must have the same length.")
  }

  if (nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("'Sigma' has incompatible dimensions.")
  }

  if (nrow(lambda) != p) {
    stop("'lambda' has incompatible dimensions.")
  }

  Sigma_inv <- solve_sym_pd(Sigma, epsilon = epsilon)
  SigmaM <- make_posdef(Sigma + lambda %*% t(lambda), epsilon = epsilon)

  centered <- matrix(y - mu, ncol = 1)
  aux1 <- t(lambda) %*% Sigma_inv %*% centered
  aux2 <- 1 + t(lambda) %*% Sigma_inv %*% lambda
  arg <- as.numeric(aux1 / sqrt(aux2))

  2 * mvtnorm::dmvnorm(x = y, mean = mu, sigma = SigmaM) * stats::pnorm(arg)
}

#' Matrix-variate skew-normal log-likelihood
#'
#' Computes the total log-likelihood for complete matrix-variate skew-normal
#' data. Each matrix observation is vectorized and evaluated through the
#' corresponding multivariate skew-normal density.
#'
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param muM Location matrix of dimension \eqn{p \times q}.
#' @param AM Skewness matrix of dimension \eqn{p \times q}.
#' @param SigmaM Row covariance matrix of dimension \eqn{p \times p}.
#' @param PsiM Column covariance matrix of dimension \eqn{q \times q}.
#' @param epsilon Numerical tolerance for covariance regularization.
#'
#' @return Numeric scalar containing the total log-likelihood.
#'
#' @export
loglik_mvsn <- function(dados, muM, AM, SigmaM, PsiM, epsilon = 1e-8) {

  if (length(dim(dados)) != 3L) {
    stop("'dados' must be a 3D array.")
  }

  n <- dim(dados)[3]
  mu_vec <- vec_col(muM)
  lambda_vec <- vec_col(AM)
  Sigma_full <- make_posdef(kronecker(PsiM, SigmaM), epsilon = epsilon)

  suma1 <- 0

  for (j in seq_len(n)) {
    dens <- dmvsn(
      y = vec_col(dados[, , j]),
      mu = mu_vec,
      Sigma = Sigma_full,
      lambda = lambda_vec,
      epsilon = epsilon
    )

    if (!is.finite(dens) || dens <= 0) {
      dens <- .Machine$double.xmin
    }

    suma1 <- suma1 + log(dens)
  }

  suma1
}

#' Censored matrix-variate skew-normal log-likelihood
#'
#' Computes the observed-data log-likelihood for censored matrix-variate
#' skew-normal data. Fully observed observations are evaluated by the
#' extended skew-normal density, fully censored observations by the corresponding
#' probability, and partially censored observations by conditional decomposition.
#'
#' @param cc Censoring indicator array. Entries equal to `1` indicate censored or
#' missing values.
#' @param LS Array of upper censoring limits.
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param muM Location matrix of dimension \eqn{p \times q}.
#' @param SigmaM Row covariance matrix of dimension \eqn{p \times p}.
#' @param PsiM Column covariance matrix of dimension \eqn{q \times q}.
#' @param lambdaM Skewness matrix of dimension \eqn{p \times q}.
#' @param epsilon Numerical tolerance for covariance regularization.
#'
#' @return Numeric scalar containing the observed-data log-likelihood.
#'
#' @export
loglik_mvsnc <- function(cc, LS, dados, muM, SigmaM, PsiM, lambdaM, epsilon = 1e-8) {
  if (length(dim(dados)) != 3L) {
    stop("'dados' must be a 3D array.")
  }

  if (!identical(dim(dados), dim(cc)) || !identical(dim(dados), dim(LS))) {
    stop("'dados', 'cc' and 'LS' must have the same dimensions.")
  }

  if (!all(cc %in% c(0, 1))) {
    stop("'cc' must contain only 0 and 1.")
  }

  p <- dim(dados)[1]
  q <- dim(dados)[2]
  m <- dim(dados)[3]

  ver <- numeric(m)

  Vari <- make_posdef(kronecker(PsiM, SigmaM), epsilon = epsilon)
  lambdaMaux <- matrix(vec_col(lambdaM), p * q, 1)
  Sigma <- make_posdef(Vari + lambdaMaux %*% t(lambdaMaux), epsilon = epsilon)
  Sigma_inv <- solve_sym_pd(Sigma, epsilon = epsilon)

  denom <- as.numeric(1 - t(lambdaMaux) %*% Sigma_inv %*% lambdaMaux)
  denom <- max(denom, epsilon)

  sqrtSigma <- sqrtm(Sigma, epsilon = epsilon)
  lambda <- solve(sqrtSigma, lambdaMaux) / sqrt(denom)
  varphi <- solve(sqrtSigma, lambda)
  mu <- vec_col(muM)

  for (j in seq_len(m)) {
    cc1 <- vec_col(cc[, , j])
    LS1 <- vec_col(LS[, , j])
    y1 <- vec_col(dados[, , j])

    miss <- which(cc1 == 1)
    obs <- which(cc1 == 0)

    if (length(miss) == 0L) {
      val <- MomTrunc::dmvESN(x = y1,
                              mu = mu,
                              Sigma = Sigma,
                              lambda = lambda,
                              tau = 0)

      if (!is.finite(val) || val <= 0) {
        val <- .Machine$double.xmin
      }

      ver[j] <- val
      next
    }

    if (length(miss) == p * q) {
      val <- MomTrunc::pmvESN(
        lower = y1,
        upper = LS1,
        mu = mu,
        Sigma = Sigma,
        lambda = lambda,
        tau = 0
      )

      if (!is.finite(val) || val <= 0) {
        val <- .Machine$double.xmin
      }

      ver[j] <- val
      next
    }

    Sigma_oo <- make_posdef(Sigma[obs, obs, drop = FALSE], epsilon = epsilon)
    Sigma_om <- Sigma[obs, miss, drop = FALSE]
    Sigma_mo <- Sigma[miss, obs, drop = FALSE]
    Sigma_mm <- Sigma[miss, miss, drop = FALSE]
    Sigma_oo_inv <- solve_sym_pd(Sigma_oo, epsilon = epsilon)

    muc <- mu[miss] + Sigma_mo %*% Sigma_oo_inv %*% (y1[obs] - mu[obs])
    Sc <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om
    Sc <- make_posdef(Sc, epsilon = epsilon)

    varphi_obs <- varphi[obs, , drop = FALSE]
    varphi_miss <- varphi[miss, , drop = FALSE]

    tilvarphi_o <- varphi_obs + Sigma_oo_inv %*% Sigma_om %*% varphi_miss
    c_oc <- as.numeric(1 / sqrt(1 + t(varphi_miss) %*% Sc %*% varphi_miss))
    tau_co <- as.numeric(t(tilvarphi_o) %*% (y1[obs] - mu[obs]))

    dens_obs <- MomTrunc::dmvESN(
      x = y1[obs],
      mu = mu[obs],
      Sigma = Sigma_oo,
      lambda = c_oc * sqrtm(Sigma_oo, epsilon = epsilon) %*% tilvarphi_o,
      tau = 0
    )

    prob_miss <- MomTrunc::pmvESN(
      lower = y1[miss],
      upper = LS1[miss],
      mu = as.vector(muc),
      Sigma = Sc,
      lambda = sqrtm(Sc, epsilon = epsilon) %*% varphi_miss,
      tau = tau_co
    )

    val <- dens_obs * prob_miss

    if (!is.finite(val) || val <= 0) {
      val <- .Machine$double.xmin
    }

    ver[j] <- val
  }

  sum(log(ver))
}

#' ECM estimation for the matrix-variate skew-normal model
#'
#' Fits a complete-data matrix-variate skew-normal model by an ECM algorithm.
#'
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#' @param epsilon Numerical tolerance for covariance regularization.
#'
#' @return An object of class `"MVSN"` containing:
#' \describe{
#'   \item{mu}{Estimated location matrix.}
#'   \item{A}{Estimated skewness matrix.}
#'   \item{Sigma}{Estimated row covariance matrix.}
#'   \item{Psi}{Estimated column covariance matrix.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iter}{Number of iterations.}
#'   \item{converged}{Logical value indicating convergence.}
#' }
#'
#' @keywords internal
mvsn_ecm <- function(dados, precision = 1e-8, MaxIter = 50, epsilon = 1e-8) {
  if (length(dim(dados)) != 3L) {
    stop("'dados' deve ser um array 3D com dimensoes p x q x n.")
  }

  if (anyNA(dados)) {
    stop("'dados' contem NA. Trate os dados antes de rodar o ECM.")
  }

  p <- dim(dados)[1]
  q <- dim(dados)[2]
  n <- dim(dados)[3]

  loglik <- numeric(MaxIter)
  criterio <- Inf
  count <- 0L

  mu <- apply(dados, c(1, 2), mean)
  A <- apply(dados, c(1, 2), mean)
  Sigma <- diag(p)
  Psi <- diag(q)
  Vari <- kronecker(Psi, Sigma)

  while (criterio > precision && count < MaxIter) {
    count <- count + 1L

    suma00 <- 0
    suma0 <- matrix(0, p, q)
    suma1 <- matrix(0, p, q)
    suma2 <- matrix(0, q, q)
    suma3 <- matrix(0, p, p)

    mu1 <- vec_col(mu)
    A1 <- vec_col(A)
    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])

      Mtij2 <- as.numeric(1 / (1 + t(A1) %*% Vari_inv %*% A1))
      Mtij2 <- max(Mtij2, epsilon)
      Mtij <- sqrt(Mtij2)

      mutij <- as.numeric(Mtij2 * t(A1) %*% Vari_inv %*% (y1 - mu1))
      dj <- mutij / Mtij

      E <- invmills(dj)
      t1 <- as.numeric(mutij + Mtij * E)
      t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

      suma00 <- suma00 + t2
      suma0 <- suma0 + (matrix(y1, p, q) - t1 * A)
      suma1 <- suma1 + t1 * (matrix(y1, p, q) - mu)
    }

    mu <- suma0 / n
    mu1 <- vec_col(mu)

    if (abs(suma00) <= epsilon) {
      stop("Atualizacao de 'A' falhou: denominador numericamente nulo.")
    }

    A <- suma1 / as.numeric(suma00)
    A1 <- vec_col(A)

    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])

      Mtij2 <- as.numeric(1 / (1 + t(A1) %*% Vari_inv %*% A1))
      Mtij2 <- max(Mtij2, epsilon)
      Mtij <- sqrt(Mtij2)

      mutij <- as.numeric(Mtij2 * t(A1) %*% Vari_inv %*% (y1 - mu1))
      dj <- mutij / Mtij

      E <- invmills(dj)
      t1 <- as.numeric(mutij + Mtij * E)
      t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

      resid <- matrix(y1 - mu1 - t1 * A1, ncol = 1)
      omega1 <- tcrossprod(resid) + (t2 - t1^2) * (A1 %*% t(A1))
      omega <- make_posdef(as.matrix(Matrix::nearPD((omega1 + t(omega1)) / 2)$mat), epsilon = epsilon)

      L1 <- t(chol(omega))
      aux <- somaL3(L1, Sigma, Psi, epsilon = epsilon)
      suma2 <- suma2 + aux$sPsi
    }

    Psi <- normalize_cov_constraint(suma2, epsilon = epsilon)
    Vari <- kronecker(Psi, Sigma)
    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])

      Mtij2 <- as.numeric(1 / (1 + t(A1) %*% Vari_inv %*% A1))
      Mtij2 <- max(Mtij2, epsilon)
      Mtij <- sqrt(Mtij2)

      mutij <- as.numeric(Mtij2 * t(A1) %*% Vari_inv %*% (y1 - mu1))
      dj <- mutij / Mtij

      E <- invmills(dj)
      t1 <- as.numeric(mutij + Mtij * E)
      t2 <- as.numeric(mutij^2 + Mtij2 + Mtij * mutij * E)

      resid <- matrix(y1 - mu1 - t1 * A1, ncol = 1)
      omega1 <- tcrossprod(resid) + (t2 - t1^2) * (A1 %*% t(A1))
      omega <- make_posdef(as.matrix(Matrix::nearPD((omega1 + t(omega1)) / 2)$mat), epsilon = epsilon)

      L1 <- t(chol(omega))
      aux <- somaL3(L1, Sigma, Psi, epsilon = epsilon)
      suma3 <- suma3 + aux$sSigma
    }

    Sigma <- make_posdef(suma3 / (q * n), epsilon = epsilon)
    Vari <- kronecker(Psi, Sigma)

    loglik[count] <- loglik_mvsn(dados, mu, A, Sigma, Psi, epsilon = epsilon)
    criterio <- compute_ecm_criterion(loglik, count)
  }

  loglik <- loglik[seq_len(count)]

  if (count == MaxIter && criterio > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(
    mu = mu,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    loglik = loglik[count],
    BIC = BIC,
    iter = count,
    converged = (criterio <= precision)
  )

  class(obj.out) <- "MVSN.ECM"
  obj.out
}

#' ECM estimation for the censored matrix-variate skew-normal model
#'
#' Fits a censored matrix-variate skew-normal model by an ECM algorithm. The
#' procedure handles fully observed, fully censored, and partially censored
#' matrix observations using conditional moments from skew-normal and extended
#' skew-normal truncated distributions.
#'
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param cc Censoring indicator array with the same dimensions as `dados`.
#' Entries equal to `1` indicate censored or missing values.
#' @param LS Array of upper censoring limits with the same dimensions as `dados`.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#' @param epsilon Numerical tolerance for covariance regularization.
#'
#' @return An object of class `"MVSNC.ECM"` containing:
#' \describe{
#'   \item{mu}{Estimated location matrix.}
#'   \item{Sigma}{Estimated row covariance matrix.}
#'   \item{Psi}{Estimated column covariance matrix.}
#'   \item{A}{Estimated skewness matrix.}
#'   \item{loglik}{Final observed-data log-likelihood.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iter}{Number of iterations.}
#'   \item{converged}{Logical value indicating convergence.}
#' }
#'
#' @keywords internal
mvsnc_ecm <- function(dados, cc, LS, precision = 1e-6, MaxIter = 50, epsilon = 1e-8) {

  validate_mvnc_input(dados, cc, LS)

  p <- dim(dados)[1]
  q <- dim(dados)[2]
  n <- dim(dados)[3]

  loglik <- numeric(MaxIter)
  criterio <- Inf
  count <- 0L

  dados2 <- dados
  dados2[cc == 1] <- get_observed_fill_value(dados, cc)

  muM <- apply(dados2, c(1, 2), mean)
  DeltaM <- apply(dados2, c(1, 2), mean)
  PsiM <- diag(q)
  SigmaM <- diag(p)

  Gamma <- kronecker(PsiM, SigmaM)
  Delta <- matrix(vec_col(DeltaM), p * q, 1)
  Sigma <- make_posdef(Gamma + Delta %*% t(Delta), epsilon = epsilon)

  sqrtSigma <- sqrtm(Sigma, epsilon = epsilon)
  Sigma_inv <- solve_sym_pd(Sigma, epsilon = epsilon)

  denom_shape <- as.numeric(1 - t(Delta) %*% Sigma_inv %*% Delta)
  denom_shape <- max(denom_shape, epsilon)
  shape <- solve(sqrtSigma, Delta) / sqrt(denom_shape)

  mu <- vec_col(muM)

  while (criterio > precision && count < MaxIter) {
    count <- count + 1L

    sumaw2 <- 0
    suma0 <- matrix(0, p, q)
    suma1 <- matrix(0, p, q)
    suma2 <- matrix(0, q, q)
    suma3 <- matrix(0, p, p)

    SSj <- sqrtm(Sigma, epsilon = epsilon)
    iSSj <- solve(SSj, diag(nrow(SSj)))
    shape_norm2 <- sum(shape^2)
    M2 <- 1 / (1 + shape_norm2)
    M2 <- max(M2, epsilon)
    varphi <- iSSj %*% shape

    Gamma_inv <- solve_sym_pd(Gamma, epsilon = epsilon)

    for (j in seq_len(n)) {
      cc1 <- vec_col(cc[, , j])
      LS1 <- vec_col(LS[, , j])
      y1 <- vec_col(dados[, , j])

      if (sum(cc1) == 0L) {
        aux1 <- as.numeric(t(varphi) %*% (y1 - mu))
        aux2 <- as.numeric(M2 * t(Delta) %*% Gamma_inv %*% (y1 - mu))

        T1.y <- aux2 + sqrt(M2) * invmills(aux1)
        T2.y <- aux2^2 + sqrt(M2) * aux2 * invmills(aux1) + M2

        y.hat <- matrix(y1, p * q, 1)
        yy.hat <- tcrossprod(y1)
        t1.hat <- T1.y
        t2.hat <- T2.y
        ty.hat <- T1.y * y.hat
      } else if (sum(cc1) == p * q) {
        eta <- sqrt(2 / pi) / sqrt(1 + shape_norm2)

        Sigma_use <- make_posdef((Sigma + t(Sigma)) / 2, epsilon = epsilon)

        aux3 <- MomTrunc::meanvarTMD(
          lower = c(y1),
          upper = c(LS1),
          mu = mu,
          Sigma = Sigma_use,
          lambda = shape,
          dist = "SN"
        )

        y.hat <- aux3$mean
        yy.hat <- aux3$EYY

        w0.hat <- MomTrunc::onlymeanTMD(
          lower = c(y1),
          upper = c(LS1),
          mu = c(mu),
          Sigma = Gamma,
          dist = "normal"
        )

        Ltemp <- prob_opt(
          lower = c(y1),
          upper = c(LS1),
          mean = c(mu),
          sigma = Gamma,
          nu = NULL,
          uselog2 = TRUE
        )

        LLtemp <- MomTrunc::pmvSN(
          lower = c(y1),
          upper = c(LS1),
          mu = mu,
          Sigma = Sigma_use,
          lambda = shape,
          log2 = TRUE
        )

        gamma <- eta * 2^(Ltemp - LLtemp)

        val <- as.numeric(t(shape) %*% iSSj %*% (y.hat - mu))
        gamma.ap <- invmills(val)

        bad_gamma <- is.na(gamma) || !is.finite(gamma) || gamma <= 0 || abs(gamma) > 100 * abs(gamma.ap)
        gamma <- if (bad_gamma) gamma.ap else gamma

        t1.hat <- as.numeric(M2 * t(Delta) %*% Gamma_inv %*% (y.hat - mu) + sqrt(M2) * gamma)

        t2.hat <- as.numeric(
          M2^2 * t(Delta) %*% Gamma_inv %*%
            (yy.hat - 2 * y.hat %*% t(mu) + mu %*% t(mu)) %*%
            Gamma_inv %*% Delta +
            M2 +
            gamma * M2^(3 / 2) * t(Delta) %*% Gamma_inv %*% (w0.hat - mu)
        )

        ty.hat <- M2 * (yy.hat - y.hat %*% t(mu)) %*% Gamma_inv %*% Delta +
          gamma * sqrt(M2) * w0.hat
      } else {
        obs <- which(cc1 == 0)
        miss <- which(cc1 == 1)

        Sigma_oo <- make_posdef(Sigma[obs, obs, drop = FALSE], epsilon = epsilon)
        Sigma_om <- Sigma[obs, miss, drop = FALSE]
        Sigma_mo <- Sigma[miss, obs, drop = FALSE]
        Sigma_mm <- Sigma[miss, miss, drop = FALSE]
        Sigma_oo_inv <- solve_sym_pd(Sigma_oo, epsilon = epsilon)

        muc <- mu[miss] + Sigma_mo %*% Sigma_oo_inv %*% (y1[obs] - mu[obs])
        Sc <- Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% Sigma_om
        Sc <- make_posdef((Sc + t(Sc)) / 2, epsilon = epsilon)

        tilvarphi.o <- varphi[obs, , drop = FALSE] +
          Sigma_oo_inv %*% Sigma_om %*% varphi[miss, , drop = FALSE]

        tau.co <- as.numeric(t(tilvarphi.o) %*% (y1[obs] - mu[obs]))
        lambda.co <- sqrtm(Sc, epsilon = epsilon) %*% varphi[miss, , drop = FALSE]

        tautil.cc <- tau.co / sqrt(1 + sum(lambda.co^2))

        SS.cc <- sqrtm(Sc, epsilon = epsilon)
        iSS.cc <- solve(SS.cc, diag(nrow(SS.cc)))
        varphi.cc <- iSS.cc %*% lambda.co

        Delta.cc <- SS.cc %*% lambda.co / sqrt(1 + sum(lambda.co^2))
        Gamma.cc <- make_posdef(Sc - Delta.cc %*% t(Delta.cc), epsilon = epsilon)
        mub.cc <- tautil.cc * Delta.cc
        eta.cc <- invmills(tau.co, 0, sqrt(1 + sum(lambda.co^2)))

        aux3 <- MomTrunc::meanvarTMD(
          lower = c(y1[miss]),
          upper = c(LS1[miss]),
          mu = as.vector(muc),
          Sigma = Sc,
          lambda = lambda.co,
          tau = tau.co,
          dist = "ESN"
        )

        w1.hat <- aux3$mean
        w2.hat <- aux3$EYY

        Ltemp <- prob_opt(
          lower = c(y1[miss]),
          upper = c(LS1[miss]),
          mean = c(muc - mub.cc),
          sigma = Gamma.cc,
          nu = NULL,
          uselog2 = TRUE
        )

        LLtemp <- MomTrunc::pmvESN(
          lower = y1[miss],
          upper = LS1[miss],
          mu = as.vector(muc),
          Sigma = Sc,
          lambda = lambda.co,
          tau = tau.co,
          log2 = TRUE
        )

        gamma.cc <- eta.cc * 2^(Ltemp - LLtemp)

        val <- as.numeric(tau.co + t(lambda.co) %*% iSS.cc %*% (w1.hat - muc))
        gamma.cc.ap <- invmills(val)

        bad_gamma <- is.na(gamma.cc) || !is.finite(gamma.cc) || gamma.cc <= 0 ||
          abs(gamma.cc) > 100 * abs(gamma.cc.ap)
        gamma.cc <- if (bad_gamma) gamma.cc.ap else gamma.cc

        w0.hat <- MomTrunc::onlymeanTMD(
          lower = y1[miss],
          upper = LS1[miss],
          mu = c(muc - mub.cc),
          Sigma = Gamma.cc,
          dist = "normal"
        )

        y0.hat <- matrix(y1, p * q, 1)
        y0.hat[miss] <- w0.hat

        aux4 <- gamma.cc * M2^(3 / 2) * t(Delta) %*% Gamma_inv %*% (y0.hat - mu)
        aux5 <- gamma.cc * sqrt(M2) * y0.hat

        y.hat <- matrix(y1, p * q, 1)
        y.hat[miss] <- w1.hat

        yy.hat <- tcrossprod(y1)
        yy.hat[obs, miss] <- y1[obs] %*% t(w1.hat)
        yy.hat[miss, obs] <- w1.hat %*% t(y1[obs])
        yy.hat[miss, miss] <- w2.hat

        t1.hat <- as.numeric(M2 * t(Delta) %*% Gamma_inv %*% (y.hat - mu) + sqrt(M2) * gamma.cc)

        t2.hat <- as.numeric(
          M2^2 * t(Delta) %*% Gamma_inv %*%
            (yy.hat - 2 * y.hat %*% t(mu) + mu %*% t(mu)) %*%
            Gamma_inv %*% Delta +
            M2 + aux4
        )

        ty.hat <- M2 * (yy.hat - y.hat %*% t(mu)) %*% Gamma_inv %*% Delta + aux5
      }

      suma0 <- suma0 + matrix(y.hat, p, q) - t1.hat * DeltaM
      suma1 <- suma1 + matrix(ty.hat, p, q) - t1.hat * muM
      sumaw2 <- sumaw2 + t2.hat

      omega1 <- yy.hat -
        y.hat %*% t(mu) -
        ty.hat %*% t(Delta) -
        mu %*% t(y.hat) +
        mu %*% t(mu) +
        t1.hat * (mu %*% t(Delta)) -
        Delta %*% t(ty.hat) +
        t1.hat * (Delta %*% t(mu)) +
        t2.hat * (Delta %*% t(Delta))

      omega <- as.matrix(Matrix::nearPD((omega1 + t(omega1)) / 2)$mat)
      omega <- make_posdef(omega, epsilon = epsilon)

      L1 <- t(chol(omega))
      aux <- somaL3(L1, SigmaM, PsiM, epsilon = epsilon)

      suma2 <- suma2 + aux$sPsi
      suma3 <- suma3 + aux$sSigma
    }

    muM <- suma0 / n

    if (abs(sumaw2) <= epsilon) {
      stop("Atualizacao de 'DeltaM' falhou: denominador numericamente nulo.")
    }

    DeltaM <- suma1 / sumaw2
    SigmaM <- make_posdef(suma3 / (q * n), epsilon = epsilon)
    PsiM <- normalize_cov_constraint(suma2, epsilon = epsilon)

    Gamma <- kronecker(PsiM, SigmaM)
    Delta <- matrix(vec_col(DeltaM), p * q, 1)
    Sigma <- make_posdef(Gamma + Delta %*% t(Delta), epsilon = epsilon)

    sqrtSigma <- sqrtm(Sigma, epsilon = epsilon)
    Sigma_inv <- solve_sym_pd(Sigma, epsilon = epsilon)

    denom_shape <- as.numeric(1 - t(Delta) %*% Sigma_inv %*% Delta)
    denom_shape <- max(denom_shape, epsilon)
    shape <- solve(sqrtSigma, Delta) / sqrt(denom_shape)

    mu <- vec_col(muM)

    loglik[count] <- loglik_mvsnc(cc, LS, dados, muM, SigmaM, PsiM, DeltaM, epsilon = epsilon)
    criterio <- compute_ecm_criterion(loglik, count)
  }

  loglik <- loglik[seq_len(count)]

  if (count == MaxIter && criterio > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(
    mu = muM,
    Sigma = SigmaM,
    Psi = PsiM,
    A = DeltaM,
    loglik = loglik[count],
    BIC = BIC,
    iter = count,
    converged = (criterio <= precision)
  )

  class(obj.out) <- "MVSNC.ECM"
  obj.out
}


