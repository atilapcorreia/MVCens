#' Generate Matrix Variate Skew Normal Samples (Complete Data)
#'
#' Generates a sample of matrices following the matrix-variate skew-normal
#' model \eqn{X = M + A |Z| + \mathcal{V}}, where \eqn{Z \sim N(0,1)} and
#' \eqn{\mathcal{V} \sim \mathcal{MVN}_{p \times q}(0, U, V)}.
#'
#' To be more precise, this function simulates a sample of size \eqn{n}
#' using the stochastic representation:
#' \deqn{X_i = A|Z_i| + \text{rmatrixnorm}(M, U, V),}
#' where \code{rmatrixnorm()} is provided by the \pkg{LaplacesDemon} package.
#'
#' @param n Integer. Number of matrices to generate.
#' @param M Mean matrix of dimension \eqn{p \times q}.
#' @param A Skewness matrix of dimension \eqn{p \times q}.
#' @param U Row covariance matrix \eqn{p \times p}.
#' @param V Column covariance matrix \eqn{q \times q}.
#'
#' @return A 3D array of generated samples with dimensions
#'   \eqn{p \times q \times n}.
#'
#' @export
#'
#' @importFrom LaplacesDemon rmatrixnorm
#' @importFrom stats rnorm
#'
#' @examples
#' p <- 3
#' q <- 2
#' n <- 20
#' M <- matrix(0, p, q)
#' A <- matrix(1, p, q)
#' U <- diag(p)
#' V <- diag(q)
#'
#' samples <- rmvsn(n = n, M = M, A = A, U = U, V = V)
#' dim(samples)   # Should be 3 x 2 x 20

rmvsn <- function(n = n,  M = M, A = A, U = U, V = V){

  # Extract the row dimension p from the covariance matrix U (p x p)
  p <- nrow(U)

  # Extract the column dimension q from the covariance matrix V (q x q)
  q <- ncol(V)

  # Pre-allocate a 3D array to store the n generated matrices (p x q x n)
  X.or <- array(data = NA, dim = c(p, q, n))

  # Loop over the sample size n
  for (i in 1:n) {

    # Generate a single latent normal variable Z ~ N(0,1)
    Xo <- rnorm(1)

    # Draw a matrix-variate normal sample with mean M and covariances U and V
    Aux <- LaplacesDemon::rmatrixnorm(M = M, U = U, V = V)

    # Construct the skew-normal sample: A * |Z| + V, where V is the MVN matrix
    X.or[ , , i] <- A * abs(Xo) + Aux
  }

  # Return the 3D array of generated samples
  return(X.or)

}

#' ECM Algorithm for the MVN, MVSN, MVNC and MVSNC models
#'
#' This function provides a unified ECM-type estimation method for several matrix-variate distributions —
#' MVN, MVSN, MVNC, and MVSNC — allowing for data sets that include interval-censored or missing values.
#'
#' @param dist A character string specifying the distribution/model to be fitted.
#'   Options are:
#'   "MVN"   – Matrix Variate Normal,
#'   "MVNC"  – Censored Matrix Variate Normal,
#'   "MVSN"  – Matrix Variate Skew-Normal,
#'   "MVSNC" – Censored Matrix Variate Skew-Normal.
#' @param samples A 3D array of dimension p × q × n containing the observed matrix-valued responses.
#' For uncensored entries, the array stores the observed values; for censored entries, it stores the lower limits..
#' @param cc A 3D indicator array with 1 for interval-censored entries and 0 for observed values.
#' @param LS A 3D array storing upper limits for interval-censored entries, and 0 for non-censored entries.
#' @param precision Convergence tolerance.
#' @param MaxIter Maximum number of iterations.
#'
#' @return
#' A list with components:
#' \item{mu}{Estimated location matrix (p × q).}
#' \item{A}{Estimated skewness matrix (p × q). Returned only for MVSN and MVSNC models.}
#' \item{Sigma}{Estimated row covariance matrix (p × p).}
#' \item{Psi}{Estimated column covariance matrix (q × q).}
#' \item{loglik}{Final log-likelihood value achieved by the ECM algorithm.}
#' \item{BIC}{Bayesian Information Criterion.}
#' \item{iter}{Number of ECM iterations performed.}
#'
#' @export
#'
#' @importFrom matrixNormal vec dmatnorm pmatnorm
#' @importFrom MomTrunc meanvarTMD
#' @importFrom mnormt sadmvn dmnorm
#' @importFrom Matrix nearPD chol
#' @importFrom stats sd dnorm pnorm
#' @importFrom MomTrunc pmvESN pmvSN onlymeanTMD
#'
#' @examples
#'
#' ########################################################
#' # Example 1: Fit an MVNC model to simulated 2x2 matrices
#' ########################################################
#'
#' set.seed(123)
#'
#' # Dimensions: p = 2, q = 2, n = 20
#' p <- 2
#' q <- 2
#' n <- 20
#'
#' # True parameters
#' M   <- matrix(0, p, q)
#' A   <- matrix(0, p, q)
#' U   <- diag(p)
#' V   <- diag(q)
#'
#' # Simulated incomplete data
#' samples <- rmatrix_censored(n = n,
#'                             cens = 0.15,
#'                             Ind = 1,
#'                             M = M,
#'                             U = U,
#'                             V = V,
#'                             A = A,
#'                             dist = "Normal")
#' dados <- samples$X.cens
#' cc    <- samples$cc
#' LS    <- samples$LS
#'
#' # Run ECM algorithm for censored MVN model
#' out <- mv_ecm("MVNC", samples = dados, cc = cc, LS = LS, MaxIter = 20)
#'
#' # Inspect estimated parameters
#' out$mu
#' out$Sigma
#' out$Psi
#' out$loglik
#' out$BIC
#' out$iter
#'
#' ##################################################################
#' # Example 2: Fit an MVN model to simulated 2x2 matrix-variate data
#' ##################################################################
#' set.seed(123)
#'
#' # Dimensions: p = 2, q = 2, n = 20
#' p <- 2
#' q <- 2
#' n <- 20
#'
#' # True parameters
#' M_true     <- matrix(0, p, q)
#' Sigma_true <- diag(p)
#' Psi_true   <- diag(q)
#'
#' # Generate sample matrices from MVN(M_true, Sigma_true, Psi_true)
#' samples <- array(0, dim = c(p, q, n))
#' for (i in 1:n) {
#'   samples[ , , i] <- LaplacesDemon::rmatrixnorm(
#'     M = M_true,
#'     U = Sigma_true,
#'     V = Psi_true
#'   )
#' }
#'
#' # Run ECM algorithm for the Matrix-Variate Normal model
#' out <- mv_ecm("MVN", samples = samples, cc = cc, LS = LS, MaxIter = 20)
#'
#' # Inspect fitted parameters
#' out$mu
#' out$Sigma
#' out$Psi
#' out$loglik
#' out$BIC
#' out$iter
#'
#' ###################################################################
#' # Example 3: Fit an MVSN model to simulated 2x2 matrix-variate data
#' ###################################################################
#' p <- 3
#' q <- 2
#' n <- 20
#' M <- matrix(0, p, q)
#' A <- matrix(1, p, q)
#' U <- diag(p)
#' V <- diag(q)
#'
#' # Generate simple test data
#'
#' set.seed(123)
#' dados <- rmvsn(n = n, M = M, A = A, U = U, V = V)
#'
#' out <- mv_ecm("MVSN", samples = dados, cc = cc, LS = LS, MaxIter = 20)
#' out$mu
#' out$A
#' out$Sigma
#' out$Psi
#' out$loglik
#' out$BIC
#' out$iter
#'
#' str(out)
#'
#' ####################################################################
#' # Example 4: Fit an MVSNC model to simulated 2x2 matrix-variate data
#' ####################################################################
#'
#' p <- 2
#' q <- 2
#' n <- 20
#' M <- matrix(0, p, q)
#' A <- matrix(1, p, q)
#' U <- diag(p)
#' V <- diag(q)
#'
#' # Simulated incomplete data
#' samples <- rmatrix_censored(n = n, cens = 0.15, Ind = 1, M = M, U = U, V = V, A = A, dist = "SN")
#' dados <- samples$X.cens
#' cc    <- samples$cc
#' LS    <- samples$LS
#'
#' out <- mv_ecm("MVSNC", samples = dados, cc = cc, LS = LS, MaxIter = 20)
#'
#' out$mu
#' out$A
#' out$Sigma
#' out$Psi
#' out$loglik
#' out$BIC
#' out$iter
#'
#' str(out)

mv_ecm <- function(dist,
                   samples = NULL,
                   cc = NULL, LS = NULL,
                   precision = 1e-7, MaxIter = 50) {

  # ---------- MVNC ----------
  if (dist == "MVNC") {
    if (is.null(samples) || is.null(cc) || is.null(LS))
      stop("For dist = 'MVNC', you must provide samples, cc and LS.")

    return(
      mvnc_ecm(samples = samples,
               cc = cc, LS = LS,
               precision = precision,
               MaxIter = MaxIter)
    )
  }

  # ---------- MVN ----------
  if (dist == "MVN") {
    if (is.null(samples))
      stop("For dist = 'MVN', you must provide samples.")

    return(
      mvn_ecm(samples = samples,
              precision = precision,
              MaxIter = MaxIter)
    )
  }

  # ---------- MVSN ----------
  if (dist == "MVSN") {
    if (is.null(samples))
      stop("For dist = 'MVSN', you must provide dados.")

    return(
      mvsn_ecm(dados = samples,
               precision = precision,
               MaxIter = MaxIter)
    )
  }

  # ---------- MVSNC ----------
  if (dist == "MVSNC") {
    if (is.null(samples) || is.null(cc) || is.null(LS))
      stop("For dist = 'MVSNC', you must provide samples, cc and LS.")

    return(
      mvsnc_ecm(dados = samples,
                cc = cc, LS = LS,
                precision = precision,
                MaxIter = MaxIter)
    )
  }

  stop("dist must be one of: 'MVNC', 'MVN', 'MVSN', 'MVSNC'.")
}
