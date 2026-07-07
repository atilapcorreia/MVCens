#' Multivariate skew-t density
#'
#' Evaluates the density of a multivariate skew-t distribution using a
#' latent-variable parameterization.
#'
#' @param y Numeric vector of observed values.
#' @param mu Numeric location vector.
#' @param Sigma Positive definite scale matrix.
#' @param lambda Numeric skewness vector.
#' @param nu Positive scalar. Degrees of freedom.
#' @param epsilon Numerical tolerance for positive-definiteness corrections.
#' @param log Logical. If `TRUE`, returns the log-density.
#'
#' @return A numeric density value, or log-density if `log = TRUE`.
#'
#' @export
dmvst <- function(y, mu, Sigma, lambda, nu, epsilon = 1e-8, log = FALSE) {

  y      <- as.vector(y)
  mu     <- as.vector(mu)
  lambda <- as.vector(lambda)

  if (!is.matrix(Sigma)) {
    Sigma <- as.matrix(Sigma)
  }

  p <- ncol(Sigma)

  if (length(y) != p || length(mu) != p || length(lambda) != p) {
    stop("'y', 'mu' and 'lambda' must have length equal to ncol(Sigma).")
  }

  if (nrow(Sigma) != p) {
    stop("'Sigma' must be a square matrix.")
  }

  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 0) {
    stop("'nu' must be a positive finite scalar.")
  }

  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  lambda_min <- min(eig$values)

  if (lambda_min <= epsilon) {

    Sigma <- Sigma + (epsilon - lambda_min) * diag(nrow(Sigma))

  }

  Sigma_inv <- chol2inv(chol(Sigma))

  centered <- matrix(y - mu, ncol = 1)
  lambda1 <- matrix(lambda, ncol = 1)

  Gama <- Sigma + lambda1 %*% t(lambda1)
  Gama <- (Gama + t(Gama)) / 2
  eig_g <- eigen(Gama, symmetric = TRUE)
  lambda_min_g <- min(eig_g$values)

  if (lambda_min_g <= epsilon) {

    Gama <- Gama + (epsilon - lambda_min_g) * diag(nrow(Gama))

  }

  Gama_inv <- chol2inv(chol(Gama))

  aux1 <- as.numeric(t(lambda1) %*% Sigma_inv %*% centered)
  aux2 <- as.numeric(1 + t(lambda1) %*% Sigma_inv %*% lambda1)
  aux2 <- max(aux2, epsilon)

  ds <- as.numeric(t(centered) %*% Gama_inv %*% centered)
  ds <- max(ds, 0)

  aux3 <- sqrt(nu + p) * aux1 / sqrt(aux2 * (nu + ds))

  log_dens <- log(2) +
    mnormt::dmt(
      x = y,
      mean = mu,
      S = Gama,
      df = nu,
      log = TRUE
    ) +
    stats::pt(aux3, df = nu + p, log.p = TRUE)

  if (!is.finite(log_dens)) {
    log_dens <- log(.Machine$double.xmin)
  }

  if (isTRUE(log)) {
    return(as.numeric(log_dens))
  }

  as.numeric(exp(log_dens))
}

#' Matrix-variate skew-t log-likelihood
#'
#' Computes the log-likelihood for complete matrix-variate skew-t data. Each
#' matrix observation is vectorized and evaluated by `dmvst()`.
#'
#' @param nu Positive scalar. Degrees of freedom.
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param muM Location matrix of dimension \eqn{p \times q}.
#' @param AM Skewness matrix of dimension \eqn{p \times q}.
#' @param SigmaM Row scale matrix of dimension \eqn{p \times p}.
#' @param PsiM Column scale matrix of dimension \eqn{q \times q}.
#' @param epsilon Numerical tolerance for positive-definiteness corrections.
#'
#' @return Numeric scalar containing the total log-likelihood.
#'
#' @export
loglik_mvst <- function(nu, dados, muM, AM, SigmaM, PsiM, epsilon = 1e-8) {
  if (length(dim(dados)) != 3L) {
    stop("'dados' must be a 3D array.")
  }

  p <- dim(dados)[1]
  q <- dim(dados)[2]
  n <- dim(dados)[3]

  if (!all(dim(muM) == c(p, q))) {
    stop("'muM' must have dimensions p x q.")
  }

  if (!all(dim(AM) == c(p, q))) {
    stop("'AM' must have dimensions p x q.")
  }

  if (!all(dim(SigmaM) == c(p, p))) {
    stop("'SigmaM' must have dimensions p x p.")
  }

  if (!all(dim(PsiM) == c(q, q))) {
    stop("'PsiM' must have dimensions q x q.")
  }

  mu1 <- as.vector(muM)
  A1 <- as.vector(AM)

  PsiM <- (PsiM + t(PsiM)) / 2
  eig_psi <- eigen(PsiM, symmetric = TRUE)
  lambda_min_psi <- min(eig_psi$values)
  if (lambda_min_psi <= epsilon) {
    PsiM <- PsiM + (epsilon - lambda_min_psi) * diag(nrow(PsiM))
  }

  SigmaM <- (SigmaM + t(SigmaM)) / 2
  eig_sigma <- eigen(SigmaM, symmetric = TRUE)
  lambda_min_sigma <- min(eig_sigma$values)
  if (lambda_min_sigma <= epsilon) {
    SigmaM <- SigmaM + (epsilon - lambda_min_sigma) * diag(nrow(SigmaM))
  }

  Vari <- kronecker(PsiM, SigmaM)
  Vari <- (Vari + t(Vari)) / 2
  eig_vari <- eigen(Vari, symmetric = TRUE)
  lambda_min_vari <- min(eig_vari$values)
  if (lambda_min_vari <= epsilon) {
    Vari <- Vari + (epsilon - lambda_min_vari) * diag(nrow(Vari))
  }

  suma1 <- 0

  for (j in seq_len(n)) {
    suma1 <- suma1 + dmvst(
      y = as.vector(dados[, , j]),
      mu = mu1,
      Sigma = Vari,
      lambda = A1,
      nu = nu,
      epsilon = epsilon,
      log = TRUE
    )
  }

  as.numeric(suma1)
}


#' ECM estimation for the matrix-variate skew-t model
#'
#' Fits a complete-data matrix-variate skew-t model by an ECM algorithm. The
#' model is based on the stochastic representation
#'
#' \deqn{
#' X_i = \mu + U_i^{-1/2}(W_i A + V_i),
#' }
#'
#' where \eqn{W_i} is a half-normal latent variable,
#' \eqn{U_i \sim Gamma(\nu/2, \nu/2)}, and
#' \eqn{V_i \sim MN_{p \times q}(0, \Sigma, \Psi)}.
#'
#' @param dados Numeric array with dimensions \eqn{p \times q \times n}.
#' @param nu Positive scalar. Initial degrees-of-freedom value.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#' @param get.nu Logical. If `TRUE`, updates `nu` by one-dimensional likelihood
#' maximization at each iteration. If `FALSE`, keeps `nu` fixed.
#' @param epsilon Numerical tolerance for covariance regularization.
#' @param nu_bounds Numeric vector of length two. Lower and upper bounds used
#' when optimizing `nu`.
#'
#' @return An object of class `"MVST.ECM"` containing:
#' \describe{
#'   \item{mu}{Estimated location matrix.}
#'   \item{A}{Estimated skewness matrix.}
#'   \item{Sigma}{Estimated row scale matrix.}
#'   \item{Psi}{Estimated column scale matrix.}
#'   \item{nu}{Estimated or fixed degrees-of-freedom parameter.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iter}{Number of iterations performed.}
#'   \item{converged}{Logical value indicating whether convergence was reached.}
#' }
#'
#' @keywords internal
mvst_ecm <- function(dados,
                     nu = 4,
                     precision = 1e-8,
                     MaxIter = 50,
                     get.nu = TRUE,
                     epsilon = 1e-8,
                     nu_bounds = c(0.05, 150)) {
  if (length(dim(dados)) != 3L) {
    stop("'dados' deve ser um array 3D com dimensoes p x q x n.")
  }

  if (anyNA(dados)) {
    stop("'dados' contem NA. Trate os dados antes de rodar o ECM.")
  }

  if (!is.numeric(precision) || length(precision) != 1L || !is.finite(precision) || precision <= 0) {
    stop("'precision' deve ser um escalar positivo.")
  }

  if (!is.numeric(MaxIter) || length(MaxIter) != 1L || !is.finite(MaxIter) ||
      MaxIter <= 0 || MaxIter != as.integer(MaxIter)) {
    stop("'MaxIter' deve ser um inteiro positivo.")
  }

  if (!is.logical(get.nu) || length(get.nu) != 1L || is.na(get.nu)) {
    stop("'get.nu' deve ser TRUE ou FALSE.")
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("'epsilon' deve ser um escalar positivo.")
  }

  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 0) {
    stop("'nu' deve ser um escalar positivo.")
  }

  if (!is.numeric(nu_bounds) || length(nu_bounds) != 2L || any(!is.finite(nu_bounds)) ||
      nu_bounds[1] <= 0 || nu_bounds[1] >= nu_bounds[2]) {
    stop("'nu_bounds' deve conter dois limites finitos com 0 < lower < upper.")
  }

  p <- dim(dados)[1]
  q <- dim(dados)[2]
  n <- dim(dados)[3]
  p1 <- p * q

  loglik <- numeric(MaxIter)
  criterio <- Inf
  count <- 0L

  mu <- apply(dados, c(1, 2), mean)
  A <- apply(dados, c(1, 2), mean)
  Sigma <- diag(p)
  Psi <- diag(q)
  Vari <- kronecker(Psi, Sigma)

  compute_latent_stats <- function(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu) {
    diff <- y1 - mu1

    Mtij2 <- as.numeric(1 / (1 + t(A1) %*% Vari_inv %*% A1))
    Mtij2 <- max(Mtij2, epsilon)
    Mtij <- sqrt(Mtij2)

    mutij <- as.numeric(Mtij2 * t(A1) %*% Vari_inv %*% diff)
    Ass <- mutij / Mtij
    dj <- as.numeric(t(diff) %*% Gama_inv %*% diff)

    prob <- dmvst(y1, mu1, Vari, A1, nu)
    prob <- max(as.numeric(prob), .Machine$double.xmin)

    log_det_gama <- as.numeric(determinant(Gama, logarithm = TRUE)$modulus)

    log_E <-
      log(2) +
      (nu / 2) * log(nu) +
      lgamma((p1 + nu + 1) / 2) -
      lgamma(nu / 2) -
      ((p1 + 1) / 2) * log(pi) -
      0.5 * log_det_gama -
      log(prob) -
      ((p1 + nu + 1) / 2) * log(dj + nu + Ass^2)

    pt_arg <- sqrt((p1 + nu + 2) / (dj + nu)) * Ass
    pt_val <- stats::pt(pt_arg, df = p1 + nu + 2)
    pt_val <- max(as.numeric(pt_val), .Machine$double.xmin)

    log_u <-
      log(4) +
      (nu / 2) * log(nu) +
      lgamma((p1 + nu + 2) / 2) -
      lgamma(nu / 2) -
      (p1 / 2) * log(pi) -
      0.5 * log_det_gama -
      log(prob) -
      ((p1 + nu + 2) / 2) * log(dj + nu) +
      log(pt_val)

    E <- exp(log_E)
    u <- exp(log_u)
    ut1 <- as.numeric(mutij * u + Mtij * E)
    ut2 <- as.numeric(mutij^2 * u + Mtij2 + Mtij * mutij * E)

    list(
      u = as.numeric(u),
      ut1 = ut1,
      ut2 = ut2
    )
  }

  while (criterio > precision && count < MaxIter) {
    count <- count + 1L

    mu1 <- vec_col(mu)
    A1 <- vec_col(A)

    suma_u <- 0
    suma_ut2 <- 0
    suma_mu <- matrix(0, p, q)
    suma_A <- matrix(0, p, q)
    suma_Psi <- matrix(0, q, q)
    suma_Sigma <- matrix(0, p, p)

    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)
    Gama <- make_posdef(Vari + A1 %*% t(A1), epsilon = epsilon)
    Gama_inv <- solve_sym_pd(Gama, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])
      estep <- compute_latent_stats(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

      ymat <- matrix(y1, p, q)

      suma_u <- suma_u + estep$u
      suma_ut2 <- suma_ut2 + estep$ut2
      suma_mu <- suma_mu + (estep$u * ymat - estep$ut1 * A)
      suma_A <- suma_A + estep$ut1 * (ymat - mu)
    }

    if (abs(suma_u) <= epsilon) {
      stop("Atualizacao de 'mu' falhou: denominador numericamente nulo.")
    }

    mu <- suma_mu / as.numeric(suma_u)
    mu1 <- vec_col(mu)

    if (abs(suma_ut2) <= epsilon) {
      stop("Atualizacao de 'A' falhou: denominador numericamente nulo.")
    }

    A <- suma_A / as.numeric(suma_ut2)
    A1 <- vec_col(A)

    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)
    Gama <- make_posdef(Vari + A1 %*% t(A1), epsilon = epsilon)
    Gama_inv <- solve_sym_pd(Gama, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])
      estep <- compute_latent_stats(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

      diff <- matrix(y1 - mu1, ncol = 1)
      omega1 <-
        estep$u * tcrossprod(diff) -
        estep$ut1 * diff %*% t(A1) -
        estep$ut1 * A1 %*% t(diff) +
        estep$ut2 * A1 %*% t(A1)

      omega <- make_posdef(as.matrix(Matrix::nearPD((omega1 + t(omega1)) / 2)$mat), epsilon = epsilon)
      L1 <- t(chol(omega))
      aux <- somaL3(L1, Sigma, Psi, epsilon = epsilon)
      suma_Psi <- suma_Psi + aux$sPsi
    }

    Psi <- normalize_cov_constraint(suma_Psi, epsilon = epsilon)
    Vari <- kronecker(Psi, Sigma)

    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)
    Gama <- make_posdef(Vari + A1 %*% t(A1), epsilon = epsilon)
    Gama_inv <- solve_sym_pd(Gama, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- vec_col(dados[, , j])
      estep <- compute_latent_stats(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

      diff <- matrix(y1 - mu1, ncol = 1)
      omega1 <-
        estep$u * tcrossprod(diff) -
        estep$ut1 * diff %*% t(A1) -
        estep$ut1 * A1 %*% t(diff) +
        estep$ut2 * A1 %*% t(A1)

      omega <- make_posdef(as.matrix(Matrix::nearPD((omega1 + t(omega1)) / 2)$mat), epsilon = epsilon)
      L1 <- t(chol(omega))
      aux <- somaL3(L1, Sigma, Psi, epsilon = epsilon)
      suma_Sigma <- suma_Sigma + aux$sSigma
    }

    Sigma <- make_posdef(suma_Sigma / (q * n), epsilon = epsilon)
    Vari <- kronecker(Psi, Sigma)

    if (get.nu) {
      opt_nu <- stats::optimize(
        f = function(nu_candidate) {
          loglik_mvst(nu_candidate, dados, mu, A, Sigma, Psi, epsilon = epsilon)
        },
        interval = nu_bounds,
        maximum = TRUE,
        tol = 1e-5
      )

      nu <- opt_nu$maximum
      loglik[count] <- opt_nu$objective
    } else {
      loglik[count] <- loglik_mvst(nu, dados, mu, A, Sigma, Psi, epsilon = epsilon)
    }

    criterio <- compute_ecm_criterion(loglik, count)
  }

  loglik <- loglik[seq_len(count)]

  if (count == MaxIter && criterio > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1 + as.integer(get.nu)
  BIC <- -2 * loglik[count] + npar * log(n)

  obj.out <- list(
    mu = mu,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    nu = nu,
    loglik = loglik,
    BIC = BIC,
    iter = count,
    converged = (criterio <= precision)
  )

  class(obj.out) <- "MVST.ECM"
  obj.out
}













