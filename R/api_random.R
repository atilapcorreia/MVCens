#' Generate observations from a registered matrix-variate model
#'
#' Validates the common matrix-variate generation contract and dispatches to
#' the generator registered for `model`. The supported models follow the
#' stochastic constructions used in the accompanying articles: matrix-normal
#' row/column covariance for `MVN`, censoring and missingness for `MVNC`,
#' half-normal skewing for `MVSN` and `MVSNC`, a Gamma scale mixture for
#' `MVST`, row-specific half-normal effects for `MVRSN`, and row-specific
#' exponential effects for `MVREN`. All generated observations use the same
#' `p` by `q` shape as `M`; `Sigma` is the `p` by `p` row covariance or scale
#' component and `Psi` is the `q` by `q` column component. The model name is
#' case-insensitive.
#'
#' The following models are available: `MVN`, `MVNC`, `MVSN`, `MVSNC`, `MVST`,
#' `MVRSN`, `MVREN`, `MVNIG`, and `MVVG`. `MVNC`, `MVSNC`, and `MVNIG` can
#' return censored-data structures; `MVVG` uses `rate`; `MVST` uses `nu`; and
#' the high-level `MVREN` interface uses the identified convention
#' `lambda_i = 1`.
#'
#' @param model One registered model name. Matching is case-insensitive.
#' @param n Positive integer giving the number of matrix observations.
#' @param M Numeric `p` by `q` location matrix.
#' @param A Optional `p` by `q` skewness or latent-effect matrix. It is
#'   required by `MVSN`, `MVSNC`, `MVST`, `MVRSN`, `MVREN`, `MVNIG`, and
#'   `MVVG`; it is ignored by `MVN` and `MVNC`.
#' @param Sigma Positive-definite `p` by `p` row covariance or scale matrix.
#' @param Psi Positive-definite `q` by `q` column covariance or scale matrix.
#' @param cens Censoring proportion in `[0, 1]`. Required for `MVNC` and
#'   `MVSNC`; for `MVNIG`, `NULL` is interpreted as `0`.
#' @param Ind Integer censoring/missingness mechanism used by `MVNC`, `MVSNC`,
#'   and `MVNIG`: `1` for interval censoring, `2` for missing values, and `3`
#'   for a mixture of interval censoring and missingness.
#' @param nu Degrees of freedom for `MVST`. This argument is required when
#'   generating from `MVST`.
#' @param lambda Optional row-specific exponential rate vector for `MVREN`.
#'   Under the package identification convention, all supplied entries must be
#'   equal to one. `NULL` uses a vector of ones.
#' @param gamma_tilde Positive inverse-Gaussian parameter used by `MVNIG`.
#' @param rate Positive Gamma rate parameter used by `MVVG`.
#' @param return_latent Logical. For `MVRSN` and `MVREN`, if `TRUE`, return
#'   both the generated array `X` and the latent variables `W`; otherwise return
#'   only the generated array.
#'
#' @return The return type depends on the model and generation options:
#' \itemize{
#'   \item `MVN`, `MVSN`, `MVST`, and `MVVG` return a numeric array with
#'   dimensions `p` by `q` by `n`.
#'   \item `MVRSN` and `MVREN` also return that array by default; with
#'   `return_latent = TRUE`, they return a list with components `X` (the
#'   generated array) and `W` (the latent variables).
#'   \item `MVNC` and `MVSNC` return a list with `X.cens`, `cc`, and `LS`.
#'   \item `MVNIG` returns a list with `X.cens`, `cc`, `LS`, and `X.or`.
#' }
#' In particular, an ordinary `MVN` call returns the array itself, so it should
#' be used as `x`, not `x$X`.
#'
#' @details Random-number generation is delegated to the model specification;
#'   this function does not fit parameters or alter the supplied scientific
#'   model. Set the R random seed before calling it when reproducibility is
#'   required. The returned object should be treated as generated data, not as
#'   evidence that a model has been fitted or validated for a real data set.
#'
#' @examples
#' ## Common parameter matrices used throughout the examples
#' set.seed(123)
#' p <- 2
#' q <- 2
#' M <- matrix(0, p, q)
#' A <- matrix(c(0.4, 0.2, 0.1, 0.3), p, q)
#' Sigma <- matrix(c(1.0, 0.2, 0.2, 1.0), p, p)
#' Psi <- matrix(c(1.0, 0.1, 0.1, 1.0), q, q)
#'
#' ## 1. MVN: complete matrix-normal data
#' x_mvn <- mv_random(
#'   model = "MVN", n = 5, M = M,
#'   Sigma = Sigma, Psi = Psi
#' )
#' dim(x_mvn)
#' # MVN returns the array itself, not a list: use x_mvn, not x_mvn$X.
#'
#' ## 2. MVNC: matrix-normal data with the three supported mechanisms
#' x_mvnc_interval <- mv_random(
#'   model = "MVNC", n = 5, M = M,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 1
#' )
#' x_mvnc_missing <- mv_random(
#'   model = "MVNC", n = 5, M = M,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 2
#' )
#' x_mvnc_mixed <- mv_random(
#'   model = "MVNC", n = 5, M = M,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 3
#' )
#' names(x_mvnc_interval)
#'
#' ## 3. MVSN: complete matrix-variate skew-normal data
#' x_mvsn <- mv_random(
#'   model = "MVSN", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi
#' )
#' dim(x_mvsn)
#'
#' ## 4. MVSNC: skew-normal data with each censoring/missingness mechanism
#' x_mvsnc_interval <- mv_random(
#'   model = "MVSNC", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 1
#' )
#' x_mvsnc_missing <- mv_random(
#'   model = "MVSNC", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 2
#' )
#' x_mvsnc_mixed <- mv_random(
#'   model = "MVSNC", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 3
#' )
#' names(x_mvsnc_interval)
#'
#' ## 5. MVST: matrix-variate skew-t data; nu is required
#' x_mvst <- mv_random(
#'   model = "MVST", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, nu = 6
#' )
#' dim(x_mvst)
#'
#' ## 6. MVRSN: row skew-normal data
#' x_mvrsn <- mv_random(
#'   model = "MVRSN", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi
#' )
#' mvrsn_latent <- mv_random(
#'   model = "MVRSN", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, return_latent = TRUE
#' )
#' names(mvrsn_latent)
#' dim(mvrsn_latent$X)
#' dim(mvrsn_latent$W)
#'
#' ## 7. MVREN: row exponential-normal data
#' ## The identified high-level interface uses lambda_i = 1.
#' x_mvren <- mv_random(
#'   model = "MVREN", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi
#' )
#' mvren_latent <- mv_random(
#'   model = "MVREN", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, lambda = rep(1, p),
#'   return_latent = TRUE
#' )
#' names(mvren_latent)
#' dim(mvren_latent$X)
#' dim(mvren_latent$W)
#'
#' ## 8. MVNIG: complete data (cens = NULL or 0) and censored data
#' x_mvnig_complete <- mv_random(
#'   model = "MVNIG", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, gamma_tilde = 2
#' )
#' x_mvnig_interval <- mv_random(
#'   model = "MVNIG", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 1,
#'   gamma_tilde = 2
#' )
#' x_mvnig_missing <- mv_random(
#'   model = "MVNIG", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 2,
#'   gamma_tilde = 2
#' )
#' x_mvnig_mixed <- mv_random(
#'   model = "MVNIG", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, cens = 0.20, Ind = 3,
#'   gamma_tilde = 2
#' )
#' names(x_mvnig_complete)
#'
#' ## 9. MVVG: matrix-variate variance-gamma data
#' x_mvvg <- mv_random(
#'   model = "MVVG", n = 5, M = M, A = A,
#'   Sigma = Sigma, Psi = Psi, rate = 1.5
#' )
#' dim(x_mvvg)
#' @family MVCens main interface
#' @seealso [mv_fit()], [mv_monte_carlo()]
#' @export
mv_random <- function(model, n, M, A = NULL, Sigma, Psi,
                      cens = NULL, Ind = 1L, nu = NULL, lambda = NULL,
                      gamma_tilde = 2, rate = 1,
                      return_latent = FALSE) {
  spec <- get_model_spec(model)
  n <- validate_positive_integer(n, "n")

  switch(
    spec$name,
    MVN = {
      spec$validate(M = M, Sigma = Sigma, Psi = Psi, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi)
    },
    MVNC = {
      spec$validate(M = M, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind)
    },
    MVSN = {
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi)
    },
    MVSNC = {
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind)
    },
    MVST = {
      if (is.null(nu)) {
        stop("'nu' must be supplied when model = 'MVST'.", call. = FALSE)
      }
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    nu = nu, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi, nu = nu)
    },
    MVRSN = {
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    return_latent = return_latent)
    },
    MVREN = {
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    lambda = lambda, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    lambda = lambda, return_latent = return_latent)
    },
    MVNIG = {
      if (is.null(cens)) cens <- 0
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind, gamma_tilde = gamma_tilde,
                    mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    cens = cens, Ind = Ind, gamma_tilde = gamma_tilde)
    },
    MVVG = {
      spec$validate(M = M, A = A, Sigma = Sigma, Psi = Psi,
                    rate = rate, mode = "generate")
      spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                    rate = rate)
    },
    stop(sprintf("Unsupported model '%s'.", spec$name), call. = FALSE)
  )
}
