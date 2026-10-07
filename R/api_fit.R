#' Fit a registered matrix-variate model by ECM
#'
#' Validates the common fit contract and dispatches to the model-specific
#' EM-type implementation described by the corresponding model article. `X`
#' must contain `n` matrix observations stored as a `p` by `q` by `n` array.
#' Matrix-normal and skew-normal censored fits use ECM updates built around
#' truncated moments, `MVST` uses an ECME-style heavy-tailed skew model, and
#' `MVRSN`/`MVREN` use row-specific latent-variable ECM kernels. The model name
#' is case-insensitive and the scientific update equations remain in the model
#' kernel.
#'
#' Models currently fitted by this interface are `MVN`, `MVNC`, `MVSN`,
#' `MVSNC`, `MVST`, `MVRSN`, and `MVREN`. `MVNIG` and `MVVG` are generation-only
#' models in this package and fail closed with an informative error when passed
#' to `mv_fit()`.
#'
#' @param model One estimable model name. Matching is case-insensitive.
#' @param X Numeric array of dimensions `p` by `q` by `n`, containing the
#'   matrix-valued observations.
#' @param cc Optional censoring indicator array with the same dimensions as
#'   `X`; `1` marks censored/missing entries and `0` marks observed entries.
#'   Required for `MVNC` and `MVSNC`.
#' @param LS Optional array of upper censoring limits, with the same dimensions
#'   as `X`. Required by `MVNC` and `MVSNC` together with `cc`.
#' @param precision Positive finite ECM convergence tolerance.
#' @param max_iter Positive integer maximum number of ECM iterations.
#' @param samples Deprecated alias for `X`. If `X = NULL` and `samples` is
#'   supplied, `samples` is used and a deprecation warning is issued.
#' @param epsilon Positive numerical tolerance used by `MVSN`, `MVSNC`, and
#'   `MVST`. `NULL` uses the model-specific default.
#' @param nu Degrees of freedom supplied to the `MVST` fit. `NULL` uses the
#'   model-specific default.
#' @param get.nu Logical control for `MVST`; if supplied, determines whether
#'   the degrees of freedom are updated during fitting.
#' @param nu_bounds Length-two numeric vector giving the admissible interval
#'   for the `MVST` degrees-of-freedom update. `NULL` uses the model default.
#' @param normalize_Psi Logical control used by `MVRSN` to normalize `Psi` to
#'   determinant one while preserving the Kronecker covariance scale.
#' @param M_init Optional initial location matrix for `MVRSN` and `MVREN`.
#' @param A_init Optional initial latent-effect/skewness matrix for `MVRSN` and
#'   `MVREN`.
#' @param Sigma_init Optional initial row covariance matrix for `MVRSN` and
#'   `MVREN`.
#' @param Psi_init Optional initial column covariance matrix for `MVRSN` and
#'   `MVREN`.
#' @param q_policy Character control for the MVREN `Q`-matrix handling policy;
#'   one of `"warn"`, `"strict"`, or `"regularize"`. `NULL` uses the model
#'   default.
#' @param verbose Optional logical verbosity control for `MVRSN` and `MVREN`.
#'   `NULL` uses the model-specific default.
#' @param eig_floor Optional positive eigenvalue floor used by `MVRSN`.
#'   `NULL` uses the model-specific default.
#' @param monotone_tol Optional tolerance for detecting decreases in the
#'   `MVRSN` likelihood path. `NULL` uses the model-specific default.
#' @param progress_callback Optional function used internally by `MVREN` to
#'   report iteration progress. It should accept `iteration`, `max_iter`, and
#'   `criterion`.
#'
#' @return A model-specific fit object containing estimated parameters and,
#'   where implemented, the likelihood path, iteration count, convergence
#'   flag, BIC, and numerical diagnostics.
#'
#' @details ECM iterations are kept sequential because their update order is
#'   part of each model's scientific specification. The function is
#'   observational with respect to the input data: it returns a fit object and
#'   does not modify `X`, `cc`, or `LS` in place. Generation-only models fail
#'   closed instead of silently switching to another estimator. Model-specific
#'   controls are explicit arguments; controls that do not apply to the chosen
#'   model are ignored by the dispatcher.
#'
#' @examples
#' ## The fits below use very few ECM iterations so that the examples are easy
#' ## to experiment with. For real analyses, increase max_iter as needed.
#' \donttest{
#' set.seed(123)
#' p <- 2
#' q <- 2
#' n <- 30
#' M <- matrix(0, p, q)
#' A <- matrix(c(0.4, 0.2, 0.1, 0.3), p, q)
#' Sigma <- matrix(c(1.0, 0.2, 0.2, 1.0), p, p)
#' Psi <- matrix(c(1.0, 0.1, 0.1, 1.0), q, q)
#'
#' ## 1. MVN: complete matrix-normal data
#' x_mvn <- mv_random("MVN", n, M, Sigma = Sigma, Psi = Psi)
#' fit_mvn <- mv_fit(
#'   model = "MVN", X = x_mvn,
#'   precision = 1e-6, max_iter = 2
#' )
#' names(fit_mvn)
#'
#' ## Legacy compatibility: `samples` is a deprecated alias for `X`.
#' fit_mvn_legacy <- suppressWarnings(mv_fit(
#'   model = "MVN", samples = x_mvn,
#'   precision = 1e-6, max_iter = 1
#' ))
#'
#' ## 2. MVNC: fit the observed/censored array together with cc and LS.
#' d_mvnc <- mv_random(
#'   "MVNC", n, M, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 1
#' )
#' fit_mvnc <- mv_fit(
#'   model = "MVNC", X = d_mvnc$X.cens,
#'   cc = d_mvnc$cc, LS = d_mvnc$LS,
#'   precision = 1e-6, max_iter = 2
#' )
#'
#' ## The same fitting call applies to missing (Ind = 2) and mixed (Ind = 3)
#' ## samples; only the generated cc/LS structures change.
#' d_mvnc_missing <- mv_random(
#'   "MVNC", n, M, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 2
#' )
#' fit_mvnc_missing <- mv_fit(
#'   "MVNC", X = d_mvnc_missing$X.cens,
#'   cc = d_mvnc_missing$cc, LS = d_mvnc_missing$LS,
#'   max_iter = 1
#' )
#' d_mvnc_mixed <- mv_random(
#'   "MVNC", n, M, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 3
#' )
#' fit_mvnc_mixed <- mv_fit(
#'   "MVNC", X = d_mvnc_mixed$X.cens,
#'   cc = d_mvnc_mixed$cc, LS = d_mvnc_mixed$LS,
#'   max_iter = 1
#' )
#'
#' ## 3. MVSN: epsilon controls numerical truncation calculations.
#' x_mvsn <- mv_random(
#'   "MVSN", n, M, A = A, Sigma = Sigma, Psi = Psi
#' )
#' fit_mvsn <- mv_fit(
#'   model = "MVSN", X = x_mvsn,
#'   epsilon = 1e-8, precision = 1e-6, max_iter = 2
#' )
#'
#' ## 4. MVSNC: censored/missing skew-normal data.
#' d_mvsnc <- mv_random(
#'   "MVSNC", n, M, A = A, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 1
#' )
#' fit_mvsnc <- mv_fit(
#'   model = "MVSNC", X = d_mvsnc$X.cens,
#'   cc = d_mvsnc$cc, LS = d_mvsnc$LS,
#'   epsilon = 1e-8, precision = 1e-6, max_iter = 2
#' )
#' d_mvsnc_missing <- mv_random(
#'   "MVSNC", n, M, A = A, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 2
#' )
#' fit_mvsnc_missing <- mv_fit(
#'   "MVSNC", X = d_mvsnc_missing$X.cens,
#'   cc = d_mvsnc_missing$cc, LS = d_mvsnc_missing$LS,
#'   max_iter = 1
#' )
#' d_mvsnc_mixed <- mv_random(
#'   "MVSNC", n, M, A = A, Sigma = Sigma, Psi = Psi,
#'   cens = 0.15, Ind = 3
#' )
#' fit_mvsnc_mixed <- mv_fit(
#'   "MVSNC", X = d_mvsnc_mixed$X.cens,
#'   cc = d_mvsnc_mixed$cc, LS = d_mvsnc_mixed$LS,
#'   max_iter = 1
#' )
#'
#' ## 5. MVST: keep nu fixed or update it during fitting.
#' x_mvst <- mv_random(
#'   "MVST", n, M, A = A, Sigma = Sigma, Psi = Psi, nu = 6
#' )
#' fit_mvst_fixed <- mv_fit(
#'   model = "MVST", X = x_mvst, nu = 6, get.nu = FALSE,
#'   epsilon = 1e-8, precision = 1e-6, max_iter = 2
#' )
#' fit_mvst_estimated <- mv_fit(
#'   model = "MVST", X = x_mvst, nu = 6, get.nu = TRUE,
#'   nu_bounds = c(2.01, 50), epsilon = 1e-8,
#'   precision = 1e-6, max_iter = 2
#' )
#'
#' ## 6. MVRSN: default normalization, or custom initialization/controls.
#' x_mvrsn <- mv_random(
#'   "MVRSN", n, M, A = A, Sigma = Sigma, Psi = Psi
#' )
#' fit_mvrsn <- mv_fit(
#'   model = "MVRSN", X = x_mvrsn,
#'   normalize_Psi = TRUE, precision = 1e-6, max_iter = 2
#' )
#' fit_mvrsn_custom <- mv_fit(
#'   model = "MVRSN", X = x_mvrsn,
#'   normalize_Psi = FALSE,
#'   M_init = M, A_init = A,
#'   Sigma_init = Sigma, Psi_init = Psi,
#'   eig_floor = 1e-8, monotone_tol = 1e-7,
#'   verbose = FALSE, precision = 1e-6, max_iter = 2
#' )
#'
#' ## 7. MVREN: custom initialization and all Q-matrix policies.
#' x_mvren <- mv_random(
#'   "MVREN", n, M, A = A, Sigma = Sigma, Psi = Psi,
#'   lambda = rep(1, p)
#' )
#' fit_mvren_warn <- mv_fit(
#'   model = "MVREN", X = x_mvren,
#'   M_init = M, A_init = A,
#'   Sigma_init = Sigma, Psi_init = Psi,
#'   q_policy = "warn", verbose = FALSE,
#'   precision = 1e-6, max_iter = 2
#' )
#' fit_mvren_strict <- mv_fit(
#'   model = "MVREN", X = x_mvren,
#'   q_policy = "strict", precision = 1e-6, max_iter = 1
#' )
#' fit_mvren_regularized <- mv_fit(
#'   model = "MVREN", X = x_mvren,
#'   q_policy = "regularize", precision = 1e-6, max_iter = 1
#' )
#'
#' ## Optional MVREN progress callback.
#' quiet_progress <- function(iteration, max_iter, criterion) invisible(NULL)
#' fit_mvren_callback <- mv_fit(
#'   model = "MVREN", X = x_mvren,
#'   q_policy = "warn", progress_callback = quiet_progress,
#'   max_iter = 1
#' )
#'
#' }
#' @family MVCens main interface
#' @seealso [mv_random()], [mv_monte_carlo()]
#' @export
mv_fit <- function(model, X = NULL, cc = NULL, LS = NULL,
                   precision = 1e-6, max_iter = 500L,
                   samples = NULL,
                   epsilon = NULL, nu = NULL, get.nu = NULL,
                   nu_bounds = NULL, normalize_Psi = NULL,
                   M_init = NULL, A_init = NULL,
                   Sigma_init = NULL, Psi_init = NULL,
                   q_policy = NULL, verbose = NULL,
                   eig_floor = NULL, monotone_tol = NULL,
                   progress_callback = NULL) {
  if (is.null(X) && !is.null(samples)) {
    X <- samples
    warning("'samples' is deprecated; use 'X'.", call. = FALSE)
  }
  if (is.null(X)) stop("'X' must be supplied.", call. = FALSE)

  spec <- get_model_spec(model, require_fit = TRUE)

  run_ecm_model(
    spec = spec, X = X, cc = cc, LS = LS,
    precision = precision, max_iter = max_iter,
    epsilon = epsilon, nu = nu, get.nu = get.nu,
    nu_bounds = nu_bounds, normalize_Psi = normalize_Psi,
    M_init = M_init, A_init = A_init,
    Sigma_init = Sigma_init, Psi_init = Psi_init,
    q_policy = q_policy, verbose = verbose,
    eig_floor = eig_floor, monotone_tol = monotone_tol,
    progress_callback = progress_callback
  )
}
