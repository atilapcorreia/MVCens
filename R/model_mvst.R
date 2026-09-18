#' Density of the multivariate skew-\eqn{t} distribution
#'
#' Evaluates the probability density or log-density of the multivariate
#' skew-\eqn{t} kernel used by the matrix-variate skew-\eqn{t} model. The
#' distribution combines a skew-normal component with a Gamma scale mixture,
#' allowing both asymmetry and heavy tails. Matrix-variate skew-\eqn{t}
#' observations are handled by [loglik_mvst()], which vectorizes each `p` by
#' `q` slice before calling this density.
#'
#' @param y Numeric vector of observed values.
#' @param mu Numeric location vector with the same length as `y`.
#' @param Sigma Positive definite scale matrix with dimensions
#'   `length(y)` by `length(y)`.
#' @param lambda Numeric skewness vector with the same length as `y`.
#' @param nu Positive scalar. Degrees of freedom.
#' @param epsilon Numerical tolerance for positive-definiteness corrections.
#' @param log Logical. If `TRUE`, returns the log-density.
#'
#' @return A numeric scalar containing the density, or the log-density when
#'   `log = TRUE`.
#'
#' @examples
#' Y <- matrix(seq(-0.2, 0.9, length.out = 12), 3, 4)
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' dmvst(
#'   y = as.vector(Y), mu = rep(0, 12), Sigma = diag(12),
#'   lambda = as.vector(A), nu = 5
#' )
#' @family MVCens density functions
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
#' Computes the log-likelihood for complete matrix-variate skew-\eqn{t} data.
#' The MVST construction extends the skew-normal matrix model with a Gamma
#' scale mixture, so that `A` controls asymmetry while `nu` controls tail
#' weight. Each `p` by `q` observation in `X` is vectorized, using the
#' Kronecker row and column scale structure implied by `SigmaM` and `PsiM`, and
#' evaluated by [dmvst()].
#'
#' @param nu Positive scalar. Degrees of freedom.
#' @param X Numeric array with dimensions \eqn{p \times q \times n}.
#' @param muM Location matrix of dimension \eqn{p \times q}.
#' @param AM Skewness matrix of dimension \eqn{p \times q}.
#' @param SigmaM Row scale matrix of dimension \eqn{p \times p}.
#' @param PsiM Column scale matrix of dimension \eqn{q \times q}.
#' @param epsilon Numerical tolerance for positive-definiteness corrections.
#'
#' @return A finite numeric scalar containing the total log-likelihood.
#'
#' @examples
#' x <- array(seq(-0.2, 2.1, length.out = 24), dim = c(3, 4, 2))
#' A <- matrix(seq(0.1, 1.2, length.out = 12), 3, 4)
#' loglik_mvst(
#'   nu = 5, X = x, muM = matrix(0, 3, 4),
#'   AM = A, SigmaM = diag(3), PsiM = diag(4)
#' )
#' @family MVCens likelihood functions
#' @export
loglik_mvst <- function(nu, X, muM, AM, SigmaM, PsiM, epsilon = 1e-8) {
  if (length(dim(X)) != 3L) {
    stop("'X' must be a 3D array.")
  }

  p <- dim(X)[1]
  q <- dim(X)[2]
  n <- dim(X)[3]

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
      y = as.vector(X[, , j]),
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
#' Fits a complete-data matrix-variate skew-t model by an ECM algorithm.
#'
#' @param X Numeric array with dimensions \eqn{p \times q \times n}.
#' @param nu Degrees-of-freedom value greater than 2. The lower bound keeps the
#'   variance finite for covariance estimation.
#' @param precision Positive scalar. Convergence tolerance.
#' @param max_iter Positive integer. Maximum number of ECM iterations.
#' @param get.nu Logical. If `TRUE`, updates `nu` by one-dimensional likelihood
#' maximization at each iteration. If `FALSE`, keeps `nu` fixed.
#' @param epsilon Numerical tolerance for covariance regularization.
#' @param nu_bounds Numeric vector of length two. Lower and upper bounds used
#' when optimizing `nu`; the lower bound must be greater than 2.
#'
#' @return An object of class `"MVST.ECM"` containing:
#' \describe{
#'   \item{M, mu}{Estimated location matrix. `M` is the canonical name and
#'   `mu` is retained as a compatibility alias.}
#'   \item{A}{Estimated skewness matrix.}
#'   \item{Sigma}{Estimated row scale matrix.}
#'   \item{Psi}{Estimated column scale matrix.}
#'   \item{nu}{Estimated or fixed degrees-of-freedom parameter.}
#'   \item{loglik}{Final log-likelihood value.}
#'   \item{loglik_history}{Sequence of log-likelihood values.}
#'   \item{BIC}{Bayesian information criterion.}
#'   \item{iterations, iter}{Number of iterations performed. `iter` is a
#'   compatibility alias.}
#'   \item{converged}{Logical value indicating whether convergence was reached.}
#'   \item{criterion}{Final relative log-likelihood change.}
#'   \item{monotone}{Whether the likelihood history avoided material decreases.}
#'   \item{monotone_drops}{Recorded likelihood decreases beyond tolerance.}
#'   \item{normalize_Psi}{Always `TRUE`; `Psi` uses the normalized representation.}
#'   \item{npar}{Parameter count used in the BIC.}
#' }
#'
#' @keywords internal
mvst_ecm <- function(X,
                     nu = 4,
                     precision = 1e-8,
                     max_iter = 50,
                     get.nu = TRUE,
                     epsilon = 1e-8,
                     nu_bounds = c(2.01, 150)) {
  validate_mvn_input(X)

  if (!is.numeric(precision) || length(precision) != 1L || !is.finite(precision) || precision <= 0) {
    stop("'precision' deve ser um escalar positivo.")
  }

  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter <= 0 || max_iter != as.integer(max_iter)) {
    stop("'max_iter' deve ser um inteiro positivo.")
  }

  if (!is.logical(get.nu) || length(get.nu) != 1L || is.na(get.nu)) {
    stop("'get.nu' deve ser TRUE ou FALSE.")
  }

  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
    stop("'epsilon' deve ser um escalar positivo.")
  }

  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 2) {
    stop("'nu' deve ser um escalar maior que 2.")
  }

  if (!is.numeric(nu_bounds) || length(nu_bounds) != 2L || any(!is.finite(nu_bounds)) ||
      nu_bounds[1] <= 2 || nu_bounds[1] >= nu_bounds[2]) {
    stop("'nu_bounds' deve conter dois limites finitos com 2 < lower < upper.")
  }

  state <- initialize_ecm_state(X, max_iter)
  p <- state$p
  q <- state$q
  n <- state$n
  p1 <- p * q

  loglik <- state$loglik
  criterio <- state$criterion
  count <- state$iteration

  mu <- state$mu
  A <- state$A
  Sigma <- state$Sigma
  Psi <- state$Psi
  Vari <- state$Vari

  mvst_compute_latent_statistics <- function(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu) {
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

  while (criterio > precision && count < max_iter) {
    count <- count + 1L

    mu1 <- matrix_vectorize(mu)
    A1 <- matrix_vectorize(A)

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
      y1 <- matrix_vectorize(X[, , j])
      estep <- mvst_compute_latent_statistics(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

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
    mu1 <- matrix_vectorize(mu)

    if (abs(suma_ut2) <= epsilon) {
      stop("Atualizacao de 'A' falhou: denominador numericamente nulo.")
    }

    A <- suma_A / as.numeric(suma_ut2)
    A1 <- matrix_vectorize(A)

    Vari_inv <- solve_sym_pd(Vari, epsilon = epsilon)
    Gama <- make_posdef(Vari + A1 %*% t(A1), epsilon = epsilon)
    Gama_inv <- solve_sym_pd(Gama, epsilon = epsilon)

    for (j in seq_len(n)) {
      y1 <- matrix_vectorize(X[, , j])
      estep <- mvst_compute_latent_statistics(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

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
      y1 <- matrix_vectorize(X[, , j])
      estep <- mvst_compute_latent_statistics(y1, mu1, A1, Vari, Vari_inv, Gama, Gama_inv, nu)

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
          loglik_mvst(nu_candidate, X, mu, A, Sigma, Psi, epsilon = epsilon)
        },
        interval = nu_bounds,
        maximum = TRUE,
        tol = 1e-5
      )

      nu <- opt_nu$maximum
      loglik[count] <- opt_nu$objective
    } else {
      loglik[count] <- loglik_mvst(nu, X, mu, A, Sigma, Psi, epsilon = epsilon)
    }

    criterio <- compute_ecm_criterion(loglik, count)
  }

  loglik <- loglik[seq_len(count)]
  diagnostics <- ecm_output_diagnostics(loglik, criterio)

  if (count == max_iter && criterio > precision) {
    warning("The algorithm stopped after reaching the maximum number of iterations without convergence.")
  }

  npar <- 2 * (p * q) + (p * (p + 1) / 2) + (q * (q + 1) / 2) - 1 + as.integer(get.nu)
  BIC <- -2 * diagnostics$loglik + npar * log(n)

  obj.out <- list(
    M = mu,
    mu = mu,
    A = A,
    Sigma = Sigma,
    Psi = Psi,
    nu = nu,
    loglik = diagnostics$loglik,
    loglik_history = diagnostics$loglik_history,
    BIC = BIC,
    iterations = count,
    iter = count,
    converged = (criterio <= precision),
    criterion = diagnostics$criterion,
    monotone = diagnostics$monotone,
    monotone_drops = diagnostics$monotone_drops,
    normalize_Psi = TRUE,
    npar = npar
  )

  class(obj.out) <- "MVST.ECM"
  obj.out
}

#' Monte Carlo study for the matrix-variate skew-t ECM estimator.
#' @keywords internal
#' @noRd
mvst_monte_carlo <- function(sample_sizes = c(50, 100, 200, 400),
                             replications = 200,
                             truth = default_mc_truth(skew = TRUE, nu = 4),
                             precision = 1e-6, max_iter = 50,
                             seed = 123, verbose = FALSE, workers = NULL) {
  run_model_monte_carlo("MVST", sample_sizes, replications, truth,
                        precision, max_iter, seed, workers, NULL, NULL, verbose)
}

mvcens_spec_mvst <- function() {
  new_model_spec(
    name = "MVST",
    validate = function(X = NULL, M = NULL, A = NULL, Sigma = NULL,
                        Psi = NULL, nu = NULL, mode, ...) {
      if (mode == "fit") validate_mvn_input(X)
      else validate_model_parameters(M, Sigma, Psi, A, require_A = TRUE)
      if (mode == "fit" && !is.null(nu) &&
          (!is.numeric(nu) || length(nu) != 1L ||
           !is.finite(nu) || nu <= 2)) {
        stop("'nu' must be greater than 2 for MVST fitting.", call. = FALSE)
      }
      invisible(TRUE)
    },
    initialize = function(X, max_iter = 200L, ...) initialize_ecm_state(X, max_iter),
    loglik = loglik_mvst,
    generate = function(n, M, A, Sigma, Psi, nu, ...) rmvst(n, M, A, Sigma, Psi, nu),
    parameter_count = function(p, q) model_parameter_count(p, q, skew = TRUE, extra = 1L),
    fit = function(X, cc = NULL, LS = NULL, precision, max_iter, nu = 4,
                   get.nu = TRUE, nu_bounds = c(2.01, 150), epsilon = 1e-8, ...) {
      mvst_ecm(X, nu = nu, precision = precision, max_iter = max_iter,
               get.nu = get.nu, epsilon = epsilon, nu_bounds = nu_bounds)
    }
  )
}
